import enum
import math
import numpy as np
import tskit
from itertools import chain
from itertools import product
import scipy.integrate as integrate
from scipy.special import logsumexp

class Phase:
  # t: beginning time of Phase
  # n: list of coalescent rates, each could be a number, a (initial_rate, growth_rate) tuple, or any function f(t)
  # P: transition matrix at beginning of Phase
  def __init__(self, t, n, P = None):
    self.t = t
    self.t_end = math.inf
    
    if P is not None:
      assert P.shape[1] == len(n)
      self.P = P
    else:
      self.P = np.identity(len(n))
    
    with np.errstate(divide='ignore'):
      logP = np.log(self.P)
      logP[np.isneginf(logP)] = 0
      self.logP = logP
    
    self.ins_outs_mappings()
    
    self.n = []
    self.intn = []
    for n_ in n:
      if type(n_) is float:
        def func_n(x):
          return n_
        def func_intn(a, b):
          return n_ * (b - a)
      
      elif (type(n_) is tuple and n_[1] == 0):
        n_ = n_[0]
        def func_n(x):
          return n_
        def func_intn(a, b):
          return n_ * (b - a)
      
      elif type(n_) is tuple:
        n_0, r = n_
        def func_n(x):
          return n_0 * math.exp(r * x)
        def func_intn(a, b):
          return n_0 / r * (math.exp(b*r) - math.exp(a*r))
      
      elif callable(n_):
        def func_n(x):
          return n_(x)
        def func_intn(a, b):
          return integrate.quad(n_, a, b)
      
      else:
        raise Exception("not supported n type: not a number, a tuple or a function")
      
      self.n.append(func_n)
      self.intn.append(func_intn)
  
  def ins_outs_mappings(self):
    ins, outs = np.where(self.P > 0)
    self.ins2outs = {in_:set() for in_ in range(self.P.shape[0])}
    self.outs2ins = {out:set() for out in range(self.P.shape[1])}
    for in_, out in zip(ins, outs):
      self.ins2outs[in_].add(out)
      self.outs2ins[out].add(in_)
  
  # N: list of initial samples in each population
  # coals: list of sorted lists of coalescences in each population
  # t: start time
  # t_end: end time
  def logp(self, N, coals, t, t_end):
    buffer = 0
    for pop in range(len(self.n)):
      N_ = N[pop]
      ts = coals[pop]
      ts_intervals = list(zip([t]+ts, ts+[t_end]))
      for t in ts:
        buffer += math.log(self.n[pop](t))
      for a, b in ts_intervals:
        buffer -= N_ * (N_-1) / 2 * self.intn(a, b)
        N_ -= 1
