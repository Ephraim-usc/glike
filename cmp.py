import math
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
  
  Q_ = Q.copy()
  np.fill_diagonal(Q_, 0)
  
  for pop in pops:
    QQ.loc[pop+pop, pop] = ns[int(pop)]
    QQ.loc[pop+pop, pop+pop] = - ns[int(pop)] - 2 * Q_[int(pop)].sum()
    for pop2 in set(pops).difference(pop):
      QQ.loc[pop+pop, pop+pop2] = Q[int(pop), int(pop2)]
      QQ.loc[pop+pop, pop2+pop] = Q[int(pop), int(pop2)]
      QQ.loc[pop2+pop, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop] = Q[int(pop2), int(pop)]
      QQ.loc[pop+pop2, pop+pop2] = - Q_[int(pop)].sum() - Q_[int(pop2)].sum()
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
    PP.loc[pop, pop] = 1
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
    self.states = [str(i) for i in range(self.N)]
  
  def get_P(self, start, end):
    P = np.identity(self.N)
    for i, time in enumerate(self.times):
      if type(time) == list: # continuous process
        if (time[1] <= start or time[0] > end): continue
        duration = min(time[1], end) - max(time[0], start)
        P = np.matmul(P, expm(self.Qs[i] * duration))
      if type(time) != list: # mass migration
        if (time <= start or time > end): continue
        P = np.matmul(P, self.Qs[i])
    P = pd.DataFrame(P)
    P.index = P.columns = self.states
    return P
  
  def get_Q(self, start):
    P = np.zeros([self.N, self.N])
    for i, time in enumerate(self.times):
      if type(time) == list: # continuous process
        if time[1] <= start: continue
        P = self.Qs[i]
      if type(time) != list: # mass migration
        continue
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


# Coalescent likelihood of node
def CLN(scmp, pcmp, tree, node, sample_pops, ps = {}):
  pops = scmp.states
  
  children = tree.children(node)
  if len(children) == 0:
    buffer = {pop:float(sample_pops[node] == pop) for pop in pops}
    ps[node] = buffer
    return buffer
  
  if tree.time(children[0]) <= tree.time(children[1]):
    child_1, child_2 = children
  else:
    child_2, child_1 = children
  
  if child_1 not in ps:
    CLN(scmp, pcmp, tree, child_1, sample_pops, ps)
  if child_2 not in ps:
    CLN(scmp, pcmp, tree, child_2, sample_pops, ps)
  
  P = scmp.get_P(tree.time(child_1), tree.time(child_2))
  PP = pcmp.get_P(tree.time(child_2), tree.time(node))
  QQ = pcmp.get_Q(tree.time(node))
  buffer = {}
  for pop in pops:
    p = 0
    for pop_1 in pops:
      for pop_2 in pops:
        for pop_1_ in pops:
          p += PP.loc[pop_1_ + pop_2, pop + pop] * P.loc[pop_1, pop_1_] * ps[child_1][pop_1] * ps[child_2][pop_2]
    buffer[pop] = p * QQ.loc[pop + pop, pop]
  
  ps[node] = buffer
  return buffer


def log10(x):
  if x < 1e-100:
    return -math.inf
  else:
    return math.log10(x)

def log_add(x, y):
  max_ = max(x, y)
  if max_ == -math.inf:
    return -math.inf
  
  x_ = x - max_
  y_ = y - max_
  buffer = max_ + log10(math.pow(10, x_) + math.pow(10, y_))
  return buffer


def logCLN(scmp, pcmp, tree, node, sample_pops, logps = {}):
  pops = scmp.states
  
  children = tree.children(node)
  if len(children) == 0:
    buffer = {pop:log10(sample_pops[node] == pop) for pop in pops}
    logps[node] = buffer
    return buffer
  
  if tree.time(children[0]) <= tree.time(children[1]):
    child_1, child_2 = children
  else:
    child_2, child_1 = children
  
  if child_1 not in logps:
    logCLN(scmp, pcmp, tree, child_1, sample_pops, logps)
  if child_2 not in logps:
    logCLN(scmp, pcmp, tree, child_2, sample_pops, logps)
  
  P = scmp.get_P(tree.time(child_1), tree.time(child_2))
  PP = pcmp.get_P(tree.time(child_2), tree.time(node))
  QQ = pcmp.get_Q(tree.time(node))
  buffer = {}
  for pop in pops:
    logp = -math.inf
    for pop_1 in pops:
      for pop_2 in pops:
        for pop_1_ in pops:
          logp = log_add(logp, log10(PP.loc[pop_1_ + pop_2, pop + pop]) + log10(P.loc[pop_1, pop_1_]) + logps[child_1][pop_1] + logps[child_2][pop_2])
    buffer[pop] = logp + log10(QQ.loc[pop + pop, pop])
  
  logps[node] = buffer
  return buffer





demography = msprime.Demography()
demography.add_population(name="A", initial_size=1_000)
demography.add_population(name="B", initial_size=1_000)
demography.add_population(name="C", initial_size=1_000)
demography.set_migration_rate(source="A", dest="B", rate=1e-5)
demography.set_migration_rate(source="B", dest="A", rate=1e-5)
demography.add_population_split(time=200, derived=["A", "B"], ancestral="C")
trees = msprime.sim_ancestry(samples={"A": 500, "B": 500}, demography=demography, ploidy = 1)
tree = trees.first()

sample_pops = ["0"] * 500 + ["1"] * 500
nss = [[0.001, 0.001, 0.001], [0.001, 0.001, 0.001], [0.001, 0.001, 0.001]]
Qs = [np.array([[-1e-5, 1e-5, 0], [1e-5, -1e-5, 0], [0, 0, 0]]), 
      np.array([[0, 0, 1], [0, 0, 1], [0, 0, 1]]), 
      np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])]

results = []
for time in range(150, 250, 5):
  print(time)
  times = [[0, time], time, [time, 1e6]]
  scmp = CMP(Qs, nss, times)
  pcmp = scmp.s2p()
  results.append(logCLN(scmp, pcmp, tree, 1998, sample_pops, logps = {})['2'])

import matplotlib.pyplot as plt
plt.scatter(range(100, 300, 10), results)
plt.show()





demography = msprime.Demography()
demography.add_population(name="A", initial_size=1000)
demography.add_population(name="B", initial_size=1000)
demography.add_population(name="ADMIX", initial_size=1000)
demography.add_population(name="ANC", initial_size=1000)
demography.set_migration_rate(source="A", dest="B", rate=1e-5)
demography.set_migration_rate(source="B", dest="A", rate=1e-5)
demography.add_admixture(time=100, derived="ADMIX", ancestral=["A", "B"], proportions=[0.2, 0.8])
demography.add_population_split(time=500, derived=["A", "B"], ancestral="ANC")
trees = msprime.sim_ancestry(samples={"ADMIX": 1000}, demography=demography, ploidy = 1)
tree = trees.first()

sample_pops = ["2"] * 1000
nss = [[0.001, 0.001, 0.001, 0.001]] * 5
Qs = [np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]), 
      np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0.2, 0.8, 0, 0], [0, 0, 0, 1]]), 
      np.array([[-1e-5, 1e-5, 0, 0], [1e-5, -1e-5, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]), 
      np.array([[0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]]),
      np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])]

results = []
for time in range(10, 400, 10):
  print(time)
  times = [[0, time], time, [time, 500], 500, [500, 1e6]]
  scmp = CMP(Qs, nss, times)
  pcmp = scmp.s2p()
  results.append(logCLN(scmp, pcmp, tree, 1998, sample_pops, logps = {})['3'])

import matplotlib.pyplot as plt
plt.scatter(range(10, 400, 10), results)
plt.show()
