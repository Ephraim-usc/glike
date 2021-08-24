import numpy as np
import pandas as pd
from scipy.linalg import expm


#ns = [0.01, 0.02, 0.03]
#Q = np.array([[0, 0.001, 0], [0.001, 0, 0.001], [0, 0.001, 0]])

def Q2QQ(Q, ns):
  N = Q.shape[0]
  pops = [str(i) for i in range(N)]
  states = pops + [x+y for x in pops for y in pops]
  
  QQ = pd.DataFrame(np.zeros([N*N+N, N*N+N]))
  QQ.index = QQ.columns = states
  
  for pop in pops:
    QQ.loc[pop+pop, pop] = ns[int(pop)]
    QQ.loc[pop+pop, pop+pop] = - ns[int(pop)] - 2 * Q[int(pop)].sum()
    for pop2 in set(pops).difference(pop):
      QQ.loc[pop+pop, pop+pop2] = Q[int(pop), int(pop2)]
      QQ.loc[pop+pop, pop2+pop] = Q[int(pop), int(pop2)]
      QQ.loc[pop2+pop, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop2] = - Q[int(pop)].sum() - Q[int(pop2)].sum()
      for pop3 in set(pops).difference([pop, pop2]):
        QQ.loc[pop+pop2, pop+pop3] = Q[int(pop2), int(pop3)]
        QQ.loc[pop+pop2, pop3+pop2] = Q[int(pop), int(pop3)]
  
  return QQ


#P = np.array([[0,0.2,0.8],[0,1,0],[0,0,1]])

def P2PP(P):
  N = P.shape[0]
  pops = [str(i) for i in range(N)]
  states = pops + [x+y for x in pops for y in pops]
  
  PP = pd.DataFrame(np.zeros([N*N+N, N*N+N]))
  PP.index = PP.columns = states
  
  for pop in pops:
    for pop2 in pops:
      for pop3 in pops:
        for pop4 in pops:
          PP.loc[pop+pop2, pop3+pop4] = P[int(pop), int(pop3)] * P[int(pop2), int(pop4)]
  
  return PP


# Coalescent Markov Process
class CMP:
  def __init__(self, Qs, nss, times):
    self.Qs = Qs
    self.nss = nss
    self.times = times
    self.N = Qs[0].shape[0]
    self.states = [str(i) for i in range(N)]
  
  def __call__(self, start, end):
    P = np.identity(self.N)
    for i, time in enumerate(self.times):
      if type(time) == list: # continuous process
        if (time[1] <= start or time[0] > end): continue
        duration = time[1] - max(time[0], start)
        P = np.matmul(P, expm(self.Qs[i] * duration))
      if type(time) != list: # mass migration
        if (time <= start or time > end): continue
        P = np.matmul(P, self.Qs[i])
    P = pd.DataFrame(P)
    P.index = P.columns = self.states
    return P
  
  def s2p(self):
    QQs = []
    for Q, ns, time in zip(self.Qs, nss, self.times):
      if type(time) == list:
        QQs.append(Q2QQ(Q, ns).values)
      if type(time) != list:
        QQs.append(P2PP(Q).values)
    pcmp = CMP(QQs, self.nss, self.times)
    pops = [str(i) for i in range(self.N)]
    pcmp.states = pops + [x+y for x in pops for y in pops]
    return pcmp

nss = [[0.001, 0.001, 0.001], [0.001, 0.001, 0.001], [0.001, 0.001, 0.001]]
Qs = [np.array([[-0.01, 0.01, 0], [0.01, -0.01, 0], [0, 0, 0]]), 
      np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1]]), 
      np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]
times = [[0, 100], 100, [100, 1000]]
scmp = CMP(Qs, nss, times)
pcmp = scmp.s2p()


# Coalescent likelihood of genealogy
def CLG(trees, cmp, labels):



import msprime

trees = msprime.sim_ancestry(samples = 100, ploidy = 1)
tree = trees.last()

# Coalescent likelihood of node
def CLN(scmp, pcmp, tree, node, sample_pops, ps = {}):
  print(str(node))
  
  pops = scmp.states
  children = tree.children(node)
  if len(children) == 0:
    return {pop:float(sample_pops[node] == pop) for pop in pops}
  if tree.time(children[0]) <= tree.time(children[1]):
    child_1, child_2 = children
  else:
    child_2, child_1 = children
  
  if child_1 not in ps:
    ps[child_1] = CLN(scmp, pcmp, tree, child_1, sample_pops, ps)
  if child_2 not in ps:
    ps[child_2] = CLN(scmp, pcmp, tree, child_2, sample_pops, ps)
  
  P = scmp(tree.time(child_1), tree.time(child_2))
  PP = pcmp(tree.time(child_2), tree.time(node))
  buffer = {}
  for pop in pops:
    p = 0
    for pop_1 in pops:
      for pop_2 in pops:
        for pop_1_ in pops:
          p += PP.loc[pop_1_ + pop_2, pop] * P.loc[pop_1, pop_1_] * ps[child_1][pop_1] * ps[child_2][pop_2]
    buffer[pop] = p
  
  return buffer

sample_pops = ["0"] * 50 + ["1"] * 50
CLN(scmp, pcmp, tree, 190, sample_pops)

