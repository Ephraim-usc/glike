import math
import itertools
import numpy as np
import pandas as pd
import tskit
import scipy.linalg
from scipy.special import logsumexp
from scipy.linalg import expm
import msprime

import state

np.seterr(all = 'ignore')
np.set_printoptions(suppress=True)


def paired(x, log = False):
  dim1, dim2 = x.shape
  buffer = np.zeros((dim1**2, dim2**2))
  for i in range(dim1):
    for j in range(dim1):
      for k in range(dim2):
        for l in range(dim2):
          if log:
            buffer[i * dim1 + j, k * dim2 + l] = x[i, k] + x[j, l]
          else:
            buffer[i * dim1 + j, k * dim2 + l] = x[i, k] * x[j, l]
  return buffer


class Phase:
  def __init__(self, t_start, t_end):
    self.t_start = t_start
    self.t_end = t_end
    self.prev = None
    self.next = None
  
  def connect(self, phase):
    assert self.t_end == phase.t_start
    
    self.next = phase
    phase.prev = self


class DiscretePhase(Phase):
  def __init__(self, t, P):
    Phase.__init__(self, t, t)
    
    self.type = "discrete"
    self.P = P
    self.PP = paired(P)
  
  def print(self):
    print("Discrete phase at " + str(self.t_start) + " generations." )
    print(self.P)


class ContinuousPhase(Phase):
  def __init__(self, t_start, t_end, Q, n):
    Phase.__init__(self, t_start, t_end)
    
    assert Q.shape[0] == Q.shape[1]
    assert Q.shape[0] == len(n)
    
    self.type = "continuous"
    self.Q = Q
    self.n = n
    
    dim = Q.shape[0]
    QQ = np.zeros([dim ** 2, dim ** 2])
    for i in range(dim):
      for j in range(dim):
        for k in range(dim):
          for l in range(dim):
            QQ[i * dim + j, k * dim + l] = Q[i, k] * (j==l) + Q[j, l] * (i==k)
    
    for i in range(dim):
      QQ[i * dim + i, i * dim + i] -= n[i]
    
    self.QQ = QQ
  
  def get_P(self, t):
    P = expm(self.Q * t)
    return P
  
  def get_PP(self, t):
    PP = expm(self.QQ * t)
    return PP
  
  def get_transition(self, t):
    if np.abs(self.Q).max() == 0: # no migration
      dim = self.Q.shape[0]
      logP = np.full((dim, dim), -np.inf)
      np.fill_diagonal(logP, 0)
      
      logRR = np.full((dim**2, dim**2), np.nan)
      np.fill_diagonal(logRR, 0)
      for i in range(dim):
        logRR[i * dim + i, i * dim + i] = - self.n[i] * t
    
    else:
      logP = np.log(self.get_P(t))
      logPP = np.log(self.get_PP(t))
      logRR = logPP - paired(logP, log = True)
    
    transition = state.Transition(t = t, logP = logP.astype(np.double), logRR = logRR.astype(np.double))
    return transition
  
  def print(self):
    print("Continuous phase from " + str(self.t_start) + " to "  + str(self.t_end) + " generations." )
    print(self.n)
    print(self.Q)


class Demography:
  
  def __init__(self):
    self.phases = []
  
  def add_phase(self, phase):
    if len(self.phases):
      self.phases[-1].connect(phase)
    self.phases.append(phase)
  
  def at(self, t): # this algorithm always returns continuous stage
    for phase in self.phases:
      if t >= phase.t_start and t <= phase.t_end:
        return phase
    return None
  
  def get_transition(self, t1, t2):
    phase1 = self.at(t1)
    phase2 = self.at(t2)
    
    if phase1 is phase2:
      return phase1.get_transition(t2 - t1)
    
    t = t1; phase = phase1
    P = 1; PP = 1
    
    while True:
      if phase.type == "discrete":
        P = np.dot(P, phase.P)
        PP = np.dot(PP, phase.PP)
      else:
        t_ = min(phase.t_end, t2)
        P = np.dot(P, phase.get_P(t_ - t))
        PP = np.dot(PP, phase.get_PP(t_ - t))
      
      phase = phase.next
      if phase is phase2.next:
        break
      t = phase.t_start
    
    logP = np.log(P)
    logPP = np.log(PP)
    logRR = logPP - paired(logP, log = True)
    
    transition = state.Transition(t = t2 - t1, logP = logP.astype(np.double), logRR = logRR.astype(np.double))
    return transition
  
  def print(self):
    for phase in reversed(self.phases):
      phase.print()
      print(" ")


def glike(tree, demo):
  # preparing times and nodes for iteration
  times, nodes = list(zip(* sorted([(time, node) for node in tree.nodes() if \
      (time := round(tree.time(node), 5)) > 0], reverse = True)))
  time_intervals = list(zip(times, times[1:]+(0,)))
  
  # preparing roots
  time_root = times[0]
  phase_root = demo.at(time_root)
  lineages = np.array([tree.root]).astype(np.intc)
  values = np.array([range(phase_root.Q.shape[0])]).T.astype(np.intc)
  bundle_root = state.Bundle(time_root, lineages, values)
  
  # main loop
  bundle = bundle_root
  for time_interval, node in zip(time_intervals, nodes):
    parent = np.array([node]).astype(np.intc)
    children = np.array(sorted(tree.children(node), reverse = True)).astype(np.intc)
    logn = np.log(demo.at(time_interval[0]).n).astype(np.double)
    bundle = bundle.diverge(parent, children, logn)
    
    if time_interval[1] == time_interval[0]:
      continue
    
    transition = demo.get_transition(time_interval[1], time_interval[0])
    bundle = bundle.evolve(transition)
    transition.free()
  
  bundle_origin = bundle
  bundle_origin.propagate()
  
  logp = bundle_root.logp()
  bundle_root.free()
  return logp


def glike_trees(trees, demo): # trees: generator or list of trees
  logp = 0
  for tree in trees:
    logp += glike(tree, demo)
  
  return logp
