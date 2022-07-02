import math
import string
import time
import numpy as np
import pandas as pd
import scipy.linalg
from scipy.linalg import expm
from scipy.special import logsumexp
import msprime
np.seterr(divide='ignore')

def int2str(i):
    return chr(i + 65)

def str2int(s):
    return ord(s) - 48

def encode(n,N=32,D="0123456789qwertyuiopasdfghjklzxc"):
    return (encode(n//N,N)+D[n%N]).lstrip("0") if n>0 else "0"

ignore = 1e-20


class NODE:
  def __init__(self, pop1, pop2):
    self.pops = (pop1, pop2)
    self.links = [] # list that contains (u, v, pointer) tuples
  
  def propagate(self, depth):
    if depth == 0:
      self.logp = 0
    elif len(self.links) == 0:
      self.logp = - math.inf
    else:
      logps = []
      for u, _, pointer in self.links:
        if not hasattr(pointer, 'logp'):
          pointer.propagate(depth - 1)
        logps.append(pointer.logp + math.log(u))
      
      self.logp = logsumexp(logps)




class PPG: # path probability graph
  def __init__(self, lin1, lin2):
    self.lins = (lin1, lin2)
    self.links = [] # list that contains pointers to NODEs
  
  @classmethod
  def init(cls, lin1, lin2, pop1, pop2):
    ppg = PPG(lin1, lin2)
    node = NODE(pop1, pop2)
    ppg.links.append(node)
    return ppg
  
  # U: two-lineage transition matrix
  # V: two-lineage transition matrix, assuming no coalescence
  def add_layer(self, U, V):
    buffer = PPG(*self.lins)
    lul = {}
    for pointer in self.links:
      in_state = ''.join(pointer.pops)
      vector = U.loc[in_state]
      for out_state, u in vector.items():
        if out_state not in lul:
          lul[out_state] = NODE(out_state[0], out_state[1])
        v = V.loc[in_state, out_state]
        if v <= ignore: continue
        lul[out_state].links.append((u, v, pointer))
    
    for _, node in sorted(lul.items()):
      buffer.links.append(node)
    
    return buffer
  
  def logp(self, ns, depth):
    logps = []
    for pointer in self.links:
      if pointer.pops[0] != pointer.pops[1]:
        continue
      if not hasattr(pointer, 'logp'):
        pointer.propagate(depth)
      logps.append(pointer.logp + math.log(ns[pointer.pops[0]]))
    
    return logsumexp(logps)



# ab, ac, bc: PPVs
# lineage: name of the new lineage
def ppg_product(ab, ac, bc, ns, lin_new):
  if ab.lins[0] == ac.lins[0]:
    flag_ac = 0
  else:
    flag_ac = 1
  
  if ab.lins[1] == bc.lins[0]:
    flag_bc = 0
  else:
    flag_bc = 1
  
  
  dc = PPG(lin_new, ac.lins[1 ^ flag_ac])
  
  for pointer_ab in ab.links:
    if pointer_ab.pops[0] != pointer_ab.pops[1]:
      continue
    pop_a, pop_b = pointer_ab.pops
    
    for pointer_ac in ac.links:
      if pop_a != pointer_ac.pops[0 ^ flag_ac]:
        continue
      pop_c = pointer_ac.pops[1 ^ flag_ac]
      
      for pointer_bc in bc.links:
        if pop_b != pointer_bc.pops[0 ^ flag_bc]:
          continue
        if pop_c != pointer_bc.pops[1 ^ flag_bc]:
          continue
        
        n = ns[pop_a]
        node = node_product(pointer_ab, pointer_ac, pointer_bc, flag_ac, flag_bc, n)
        node.pops = (pop_a, pop_c) # merging two identical pops into pop_a
        dc.links.append(node)
  
  return dc


def node_product(ab, ac, bc, flag_ac, flag_bc, n):
  dc = NODE(''.join(ab.pops), ac.pops[1 ^ flag_ac])
  
  if len(ab.links) == 0:
    return dc
  
  for u_ab, v_ab, pointer_ab in ab.links:
    pop_a, pop_b = pointer_ab.pops
    
    for u_ac, v_ac, pointer_ac in ac.links:
      if pop_a != pointer_ac.pops[0 ^ flag_ac]:
        continue
      pop_c = pointer_ac.pops[1 ^ flag_ac]
      
      for u_bc, v_bc, pointer_bc in bc.links:
        if pop_b != pointer_bc.pops[0 ^ flag_bc]:
          continue
        if pop_c != pointer_bc.pops[1 ^ flag_bc]:
          continue
        
        node = node_product(pointer_ab, pointer_ac, pointer_bc, flag_ac, flag_bc, 1)
        u = u_ab * u_ac * u_bc / math.sqrt(v_ab * v_ac * v_bc)
        v = u_ab * math.sqrt(v_ac * v_bc / v_ab)
        
        dc.links.append((u * n, v * n, node))
  
  return dc






class Phase:
  pass


class DiscretePhase(Phase):
  def __init__(self, P, in_pops, out_pops):
    assert P.shape[0] == len(in_pops)
    assert P.shape[1] == len(out_pops)
    
    self.type = "discrete"
    self.in_pops = in_pops
    self.out_pops = out_pops
    self.in_states = [x+y for x in in_pops for y in in_pops]
    self.out_states = [x+y for x in out_pops for y in out_pops]
    self.P = pd.DataFrame(P, index = in_pops, columns = out_pops)
    
    in_dim = P.shape[0]
    out_dim = P.shape[1]
    PP = np.zeros([in_dim ** 2, out_dim ** 2])
    for i in range(in_dim):
      for j in range(in_dim):
        for k in range(out_dim):
          for l in range(out_dim):
            PP[i * in_dim + j, k * out_dim + l] = P[i, k] * P[j, l]
    
    self.PP = pd.DataFrame(PP, index = self.in_states, columns = self.out_states)
  
  def get_U(self, t):
    assert t == 0
    return self.PP
  
  def get_V(self, t):
    assert t == 0
    return self.PP


class ContinuousPhase(Phase):
  def __init__(self, Q, ns, pops):
    assert Q.shape[0] == Q.shape[1]
    assert Q.shape[0] == len(ns)
    assert Q.shape[0] == len(pops)
    
    self.type = "continuous"
    self.pops = pops
    self.states = [x+y for x in pops for y in pops]
    self.Q = pd.DataFrame(Q, index = pops, columns = pops)
    self.ns = dict(zip(pops, ns))
    
    dim = Q.shape[0]
    QQ = np.zeros([dim ** 2, dim ** 2])
    for i in range(dim):
      for j in range(dim):
        for k in range(dim):
          for l in range(dim):
            QQ[i * dim + j, k * dim + l] = Q[i, k] * (j==l) + Q[j, l] * (i==k)
    
    for i in range(dim):
      QQ[i * dim + i, i * dim + i] -= ns[i]
    
    self.QQ = pd.DataFrame(QQ, index = self.states, columns = self.states)
  
  def get_U(self, t):
    QQ = self.QQ.values
    QQexp = expm(QQ * t)
    QQexp = pd.DataFrame(QQexp, index = self.states, columns = self.states)
    return QQexp
  
  def get_V(self, t):
    Q = self.Q.values
    Qexp = expm(Q * t)
    
    dim = Q.shape[0]
    V = np.zeros([dim ** 2, dim ** 2])
    for i in range(dim):
      for j in range(dim):
        for k in range(dim):
          for l in range(dim):
            V[i * dim + j, k * dim + l] = Qexp[i, k] * Qexp[j, l]
    
    V = pd.DataFrame(V, index = self.states, columns = self.states)
    return V


class Demography:
  def __init__(self):
    self.phases = []
    self.durations = []
    self.records = {}
  
  def add_stage(self, phase, duration):
    self.phases.append(phase)
    self.durations.append(duration)
  
  def get_ns(self, t):
    assert 0 <= t
    t = round(t, 2)
    
    ns = None
    for phase, duration in zip(self.phases, self.durations):
      if t < 0:
        break
      if phase.type == "discrete":
        continue
      ns = phase.ns
      t -= duration
    
    return ns
  
  def get_UV(self, t1, t2):
    assert 0 <= t1
    assert t1 <= t2
    
    t1 = round(t1, 2)
    t2 = round(t2, 2)
    if (t1, t2) in self.records.keys():
      return self.records[(t1, t2)]
    t1_copy, t2_copy = t1, t2
    
    U = None; V = None
    for phase, duration in zip(self.phases, self.durations):
      if t2 <= 0: # when t2 == duration == 0, we assume no discrete transition
        break
      if t1 > duration: # note that t1 == duration == 0 for discrete phase
        t1 -= duration
        t2 -= duration
        continue
      
      delta = min(duration, t2) - t1
      U_ = phase.get_U(delta)
      V_ = phase.get_V(delta)
      if U is None:
        U = U_; V = V_
        index = U_.index
        columns = U_.columns
      else:
        U = np.dot(U, U_)
        V = np.dot(V, V_)
        columns = U_.columns
      
      t1 = 0
      t2 -= duration
    
    U = pd.DataFrame(U, index = index, columns = columns)
    V = pd.DataFrame(V, index = index, columns = columns)
    
    self.records[(t1_copy, t2_copy)] = (U, V)
    return U, V







def gLike(tree, demo): #currently only considering admixture population 'O'
  N = tree.num_samples()
  lins = list(tree.nodes())
  times = [tree.time(node) for node in lins]
  tmp = sorted(zip(times, lins))
  lins_all = [x[1] for x in tmp]
  lins = [None] + [x[1] for x in tmp[N:]]
  times = [0] + [x[0] for x in tmp[N:]]
  
  lins_current = list(range(N))
  ppgs = {}
  
  ppgs[0] = {}
  for i in lins_current:
    ppgs[0][i] = {}
    for j in range(i):
      ppgs[0][i][j] = PPG.init(i, j, 'O', 'O')
  
  for epoch in range(1, len(times)):
    U, V = demo.get_UV(times[epoch - 1], times[epoch])
    ns = demo.get_ns(times[epoch])
    
    ppgs[epoch] = {}
    for i in lins_current:
      ppgs[epoch][i] = {}
      for j in lins_current:
        if lins_current.index(i) <= lins_current.index(j): continue
        ppgs[epoch][i][j] = ppgs[epoch - 1][i][j].add_layer(U, V)
    
    new_lin = lins[epoch]
    children = sorted(tree.children(new_lin), key = lins_all.index)
    lins_current = [x for x in lins_current if x not in children]
    ppgs[epoch][new_lin] = {}
    for i in lins_current:
      ab = ppgs[epoch][children[1]][children[0]]
      x, y = sorted((children[1], i), key = lins_all.index)
      ac = ppgs[epoch][y][x]
      x, y = sorted((children[0], i), key = lins_all.index)
      bc = ppgs[epoch][y][x]
      
      print("epoch " + str(epoch))
      A = - time.time()
      ppgs[epoch][new_lin][i] = ppg_product(ab, ac, bc, ns = ns, lin_new = new_lin)
      A += time.time()
      print(A)
    
    lins_current.append(new_lin)
  
  #return ppgs
  a, b = sorted(tree.children(tree.root), key = lins_all.index)
  ppg = ppgs[epoch][b][a]
  logp = ppg.logp(ns, epoch)
  return ppg, logp


ppg, logp = gLike(tree, demo)






