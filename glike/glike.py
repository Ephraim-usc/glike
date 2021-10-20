import math
import numpy as np
import pandas as pd
from scipy.linalg import expm
#import pickle

def log(x):
  if x != 0:
    return math.log(x)
  else:
    return -math.inf

def mn2Q(ms, ns):
  n = len(ns)
  pops = [str(i) for i in range(n)]
  states = [x for x in pops]
  
  Q = pd.DataFrame(np.zeros([n, n]))
  Q.index = Q.columns = states
  
  for pop in pops:
    for pop2 in set(pops).difference(pop):
      Q.loc[pop, pop2] = ms[int(pop), int(pop2)]
      Q.loc[pop, pop] -= ms[int(pop), int(pop2)]
  
  return Q


def mn2QQ(ms, ns):
  n = len(ns)
  pops = [str(i) for i in range(n)]
  states = [x+y for x in pops for y in pops] + ['omega']
  
  QQ = pd.DataFrame(np.zeros([n*n+1, n*n+1]))
  QQ.index = QQ.columns = states
  
  for pop in pops:
    QQ.loc[pop+pop, "omega"] += ns[int(pop)]
    QQ.loc[pop+pop, pop+pop] = - ns[int(pop)] - 2 * ms[int(pop)].sum()
    for pop2 in set(pops).difference(pop):
      QQ.loc[pop+pop, pop+pop2] = ms[int(pop), int(pop2)]
      QQ.loc[pop+pop, pop2+pop] = ms[int(pop), int(pop2)]
      QQ.loc[pop2+pop, pop+pop] = ms[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop] = ms[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop2] = - ms[int(pop)].sum() - ms[int(pop2)].sum()
      for pop3 in set(pops).difference([pop, pop2]):
        QQ.loc[pop+pop2, pop+pop3] = ms[int(pop2), int(pop3)]
        QQ.loc[pop+pop2, pop3+pop2] = ms[int(pop), int(pop3)]
  
  return QQ


def P2PP(P):
  n = P.shape[0]
  pops = [str(i) for i in range(n)]
  states = [x+y for x in pops for y in pops] + ['omega']
  
  PP = pd.DataFrame(np.zeros([n*n+1, n*n+1]))
  PP.index = PP.columns = states
  PP.loc['omega', 'omega'] = 1
  
  for pop in pops:
    for pop2 in pops:
      for pop3 in pops:
        for pop4 in pops:
          PP.loc[pop+pop2, pop3+pop4] = P[int(pop), int(pop3)] * P[int(pop2), int(pop4)]
  
  return PP


def p_transpose(p, pops):
  p_ = p.copy()
  for pop in pops:
    for pop2 in set(pops).difference(pop):
      p_[pop + pop2] = p[pop2 + pop]
  return p_


def p_trim(p):
  p_ = p.copy()
  p_["omega"] = 0
  p_ /= p_.sum()
  return p_


def p_new(p_ab, p_ac, p_bc, pops):
  tmp = p_bc * 0
  for pop in pops:
    for pop2 in pops:
      tmp[pop + pop2] = p_bc[pop + pop]
  buffer = np.sqrt(p_ab * p_ac * tmp)
  buffer /= buffer.sum()
  return buffer


class LMP: #lineage Markov Process
  def __init__(self, times, mss, nss, Ps = None):
    self.n = len(nss[0])
    self.pops = [str(i) for i in range(self.n)]
    self.states = [x+y for x in self.pops for y in self.pops] + ["omega"]
    
    self.times = times
    self.mss = mss
    self.nss = nss
    
    if Ps is not None:
      self.Ps = Ps
    else:
      self.Ps = [np.identity(self.n)] * len(self.times)
    
    self.QQs = [mn2QQ(ms, ns) for ms, ns in zip(mss, nss)]
    self.PPs = [P2PP(P) for P in self.Ps]
  
  def get_PP(self, t1, t2):
    PP = np.identity(self.n * self.n + 1)
    for i, start in enumerate(self.times):
      end = self.times[i+1] if i+1 < len(self.times) else math.inf
      if (t2 < start or t1 >= end): continue
      if (t1 <= start):
        PP = np.matmul(PP, self.PPs[i].values)
      duration = min(t2, end) - max(t1, start)
      PP = np.matmul(PP, expm(self.QQs[i].values * duration))
    
    PP = pd.DataFrame(PP)
    PP.index = PP.columns = self.states
    return PP
  
  def get_QQ(self, t):
    QQ = None
    for i, start in enumerate(self.times):
      end = self.times[i+1] if i+1 < len(self.times) else math.inf
      if (t < start or t >= end): continue
      QQ = self.QQs[i]
    return QQ
  
  def get_P(self, t1, t2):
    P = np.identity(self.n)
    for i, start in enumerate(self.times):
      end = self.times[i+1] if i+1 < len(self.times) else math.inf
      if (t2 < start or t1 >= end): continue
      if (t1 <= start):
        P = np.matmul(P, self.Ps[i].values)
      duration = min(t2, end) - max(t1, start)
      P = np.matmul(P, expm(self.Qs[i].values * duration))
    
    P = pd.DataFrame(P)
    P.index = P.columns = self.states
    return P
  
  def get_Q(self, t):
    Q = None
    for i, start in enumerate(self.times):
      end = self.times[i+1] if i+1 < len(self.times) else math.inf
      if (t < start or t >= end): continue
      Q = self.Qs[i]
    return Q


def loglike_tree(tree, labels, lmp): # tree nodes must be sorted
  n = lmp.n
  pops = [str(i) for i in range(n)]
  states = [x+y for x in pops for y in pops] + ['omega']
  
  N = tree.num_samples()
  nodes = list(tree.nodes())
  times = [tree.time(node) for node in nodes]
  nodes = [x for _, x in sorted(zip(times, nodes))] # time sorted nodes
  
  ps = {}
  for a in range(N):
    ps[a] = {}
    for b in range(a):
      ps[a][b] = (np.array(states) == labels[a] + labels[b]).astype(float)
      ps[a][b] = pd.Series(ps[a][b], index = states)
      ps[b][a] = p_transpose(ps[a][b], pops)
  
  for a in [x for x in nodes if x not in range(N)]: # inner nodes
    ps[a] = {}
    c, d = tree.children(a)
    for b in nodes[:nodes.index(a)]:
      if tree.time(tree.parent(b)) <= tree.time(a): continue
      
      A = p_trim(np.matmul(ps[c][b], lmp.get_PP(max(tree.time(c), tree.time(b)), tree.time(a))))
      B = p_trim(np.matmul(ps[d][b], lmp.get_PP(max(tree.time(d), tree.time(b)), tree.time(a))))
      C = p_trim(np.matmul(ps[c][d], lmp.get_PP(max(tree.time(c), tree.time(d)), tree.time(a))))
      ps[a][b] = p_new(A, B, C, pops)
      ps[b][a] = p_transpose(ps[a][b], pops)
  
  Ps = {}; logP = 0
  for a in nodes:
    Ps[a] = {}
    for b in nodes[:nodes.index(a)]:
      if b not in ps[a]: continue
      
      t1 = max(tree.time(a), tree.time(b))
      t2 = min(tree.time(tree.parent(a)), tree.time(tree.parent(b)))
      p_ = np.matmul(ps[a][b], lmp.get_PP(t1, t2))
      
      if tree.parent(a) == tree.parent(b):
        Ps[a][b] = (lmp.get_QQ(t2)["omega"] * p_).values.sum()
      else:
        Ps[a][b] = 1 - p_['omega']
      
      if Ps[a][b] == 0:
        #pickle.dump(ps[a][b], open("ps.p", "wb"))
        #pickle.dump(lmp.get_PP(t1, t2), open("PP.p", "wb"))
        #pickle.dump(lmp.get_QQ(t2), open("QQ.p", "wb"))
        #pickle.dump([t1, t2], open("t1t2.p", "wb"))
      
      logP += log(Ps[a][b])
  
  return ps, Ps, logP


def loglike_trees(trees, labels, lmp, stride):
  logP = []
  for tree in trees.trees():
    if tree.index % stride != 0: continue
    print("tree " + str(tree.index))
    _, _, logP_ = loglike_tree(tree, labels, lmp)
    logP.append(logP_)
  
  return np.array(logP)




