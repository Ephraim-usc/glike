import enum
import math
import numpy as np
import tskit
from itertools import chain
from itertools import product
import scipy.integrate as integrate
from scipy.special import logsumexp
from time import time as now

def logp_n(n, a):
  if type(n) is float:
    return math.log(n)
  elif type(n) is tuple and n[1] == 0:
    n = n[0]
    return math.log(n)
  elif type(n) is tuple:
    n_0, r = n
    return math.log(n_0) + r * a
  elif callable(n):
    return math.log(n(a))
  else:
    raise Exception("not supported n type: not a number, a tuple or a function")

def logp_intn(n, a, b):
  if type(n) is float:
    return - n * (b - a)
  if type(n) is tuple and n[1] == 0:
    n = n[0]
    return - n * (b - a)
  if type(n) is tuple:
    n_0, r = n
    return - n_0 / r * (math.exp(b*r) - math.exp(a*r))
  elif callable(n):
    return - integrate.quad(n, a, b)
  else:
    raise Exception("not supported n type: not a number, a tuple or a function")


class Phase:
  # t: beginning time of Phase
  # ns: list of coalescent rates, each could be a number, a (initial_rate, growth_rate) tuple, or any function f(t)
  # P: transition matrix at beginning of Phase
  def __init__(self, t, ns, P = None):
    self.t = t
    self.t_end = math.inf
    self.ns = ns
    self.K = len(ns) # number of populations (during continuous period)
    
    if P is not None:
      assert P.shape[1] == self.K
      self.P = P
    else:
      self.P = np.identity(self.K)
    
    with np.errstate(divide='ignore'):
      logP = np.log(self.P)
      logP[np.isneginf(logP)] = 0
      self.logP = logP
    
    self.ins_outs_mappings()
  
  # build fasting mappings of non-zero migrations in P
  def ins_outs_mappings(self):
    self.ins2outs = [set() for _ in range(self.P.shape[0])]
    self.outs2ins = [set() for _ in range(self.P.shape[1])]
    ins, outs = np.where(self.P > 0)
    for in_, out in zip(ins, outs):
      self.ins2outs[in_].add(out)
      self.outs2ins[out].add(in_)
  
  # this function computes the logp of a genealogical tree during a non-migration period, given populations of lineages
  # Ns: list of initial samples in each population
  # coals: list of sorted lists of coalescences in each population
  # t: start time
  # t_end: end time
  def logp(self, Ns, coals, t, t_end):
    buffer = 0
    for pop in range(self.K):
      N = Ns[pop]
      ts = coals[pop]
      ts_intervals = list(zip([t]+ts, ts+[t_end]))
      for a in ts:
        buffer += logp_n(self.ns[pop], a-self.t)
      for a, b in ts_intervals:
        buffer += N*(N-1)/2 * logp_intn(self.ns[pop], a-self.t, b-self.t)
        N -= 1
    return buffer


class Demo:
  def __init__(self):
    self.phases = []
  
  def add_phase(self, phase):
    if len(self.phases):
      assert self.phases[-1].t < phase.t, "time error when adding phase!"
      assert self.phases[-1].K == phase.P.shape[0], "shape error when adding phase!"
      self.phases[-1].t_end = phase.t
    self.phases.append(phase)




def test_demo(t1, t2, t3, r, N_a, N_b, N_c, N_d, N_e, N_f, N_g):
  demo = Demo()
  demo.add_phase(Phase(0, np.array([N_a, N_b, N_c, N_d, N_e])))
  demo.add_phase(Phase(t1, np.array([N_a, N_c, N_d, N_e]), np.array([[1, 0, 0, 0], [r, 1-r, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])))
  demo.add_phase(Phase(t2, np.array([N_a, N_c, N_f]), np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]])))
  demo.add_phase(Phase(t3, np.array([N_g]), np.array([[1], [1], [1]])))
  return demo

demo = test_demo(30, 80, 700, 0.3, 5000, 2000, 10000, 20000, 5000, 2000, 10000)


def twoway_admixture_demo(t1, t2, r, N, N_a, N_b, N_c):
  demo = Demo()
  demo.add_phase(Phase(0, [1/N, 1/N_a, 1/N_b]))
  demo.add_phase(Phase(t1, [1/N_a, 1/N_b], P = np.array([[r, 1-r],[1, 0],[0, 1]])))
  demo.add_phase(Phase(t2, [1/N_c], P = np.array([[1], [1]])))
  return demo

def twoway_admixture_demography(t1, t2, r, N, N_a, N_b, N_c):
  import msprime
  demography = msprime.Demography()
  demography.add_population(name = "O", initial_size = N)
  demography.add_population(name = "A", initial_size = N_a)
  demography.add_population(name = "B", initial_size = N_b)
  demography.add_population(name = "C", initial_size = N_c)
  
  demography.add_admixture(time=t1, derived="O", ancestral=["A", "B"], proportions = [r, 1-r])
  demography.add_population_split(time=t2, derived=["A", "B"], ancestral="C")
  
  return demography

demo = twoway_admixture_demo(20, 500, 0.3, 2000, 10000, 30000, 5000)
demography = twoway_admixture_demography(20, 500, 0.3, 2000, 10000, 30000, 5000)






class State:
  def __init__(self):
    self.logp_evolve = math.nan # evolve log probability of this State
    self.logp = math.nan # log probability of non-coalescence and coalescence during this period
    self.children = []

class Bundle:
  # phase: the Phase this Bundle is in, the origin Bundle should have None phase
  # lins: list of lineages
  # pops: list of deterministic populations for each lineage
  def __init__(self, phase, lins, pops = None):
    self.phase = phase
    if phase is not None:
      self.t = phase.t
      self.t_end = phase.t_end
    else:
      self.t = 0
      self.t_end = 0
    
    self.lins = lins; self.N = len(lins)
    self.dict = {lin:i for i, lin in enumerate(lins)} # for fast indexing lineages
    
    self.dests = [{lin} for lin in lins] # list of sets of descendents
    self.coals = [[] for _ in lins] # list of lists of coalescent times
    if pops is not None:
      self.pops = [{pop} for pop in pops] # list of sets of possible populations
    else:
      self.pops = [{0} for _ in lins]
    
    self.states = {}
    self.parent = None
    self.child = None
  
  # backward in time through continuous demography
  # summarizing coalescent results
  def coalesce(self, t, children, parent):
    global A, B, C, D, E
    D += 1
    idx = [self.dict[child] for child in children]
    self.lins.append(parent)
    self.dests.append(set.union(*[self.dests[i] for i in idx]))
    self.pops.append(set.intersection(*[self.pops[i] for i in idx]))
    self.coals.append([coal for i in idx for coal in self.coals[i]] + [t] * (len(children) - 1))
    #print(str(parent) + " " + str(children) + " " + str(self.dests[-1]) + " " + str([self.pops[i] for i in idx]) + " " + str(self.pops[-1]))
    
    for i in sorted(idx, reverse=True):
      del self.lins[i]
      del self.dests[i]
      del self.pops[i]
      del self.coals[i]
    B -= now()
    self.dict = {lin:i for i, lin in enumerate(self.lins)}; C += 1
    B += now()
    self.N -= len(children) - 1
  
  # backward in time through mass migrations
  # creating new Bundle, computing possible populations for each lineage
  def transit(self, phase):
    bundle = Bundle(phase, self.lins.copy())
    bundle.pops = [set().union(*[phase.ins2outs[in_] for in_ in pops]) for pops in self.pops] # set().union() can handle empty pops
    
    bundle.child = self
    self.parent = bundle
    return bundle
  
  # make it the root bundle
  def root(self):
    self.t_end = max([coal for coals in self.coals for coal in coals])
    for value in product(*self.pops):
      state = State()
      self.states[value] = state
  
  # forward in time through continuous demography
  # computing coalescence and non-coalescence probabilities
  def evolve(self):
    K = self.phase.K
    for value, state in self.states.items():
      Ns = [0 for pop in range(K)]
      coals = [[] for pop in range(K)]
      for i in range(self.N):
        Ns[value[i]] += 1
        coals[value[i]].extend(self.coals[i])
      for pop in range(K):
        Ns[pop] += len(coals[pop])
        coals[pop].sort()
      state.logp_evolve = self.phase.logp(Ns, coals, self.t, self.t_end)
  
  # forward in time through mass migrations
  # creating children states and computing migration probabilities
  def migrate(self):
    outs2ins = self.phase.outs2ins
    logP = self.phase.logP
    child = self.child
    for value, state in self.states.items():
      outs = [math.nan for _ in child.lins]
      for i,pop in enumerate(value):
        for dest in self.dests[i]:
          outs[child.dict[dest]] = pop
      
      pops_child = [set.intersection(child.pops[i], outs2ins[outs[i]]) for i in range(child.N)]
      for ins in product(*pops_child):
        migrations = np.zeros(self.phase.P.shape)
        for in_, out in zip(ins, outs):
          migrations[in_, out] += 1
        logp = (logP * migrations).sum()
        
        state_child = child.states.setdefault(ins, State())
        state.children.append((logp, state_child))
  
  def origin(self):
    if self.states:
      self.logp = 0
    else:
      self.logp = -math.inf
    for state in self.states.values():
      state.logp = 0
  
  def get_logp(self):
    if len(self.states) == 0:
      self.logp = -math.inf
      return
    for state in self.states.values():
      logps = [logp + child.logp for logp, child in state.children]
      state.logp = logsumexp(logps) + state.logp_evolve
    self.logp = logsumexp([state.logp for state in self.states.values()])



def glike(tree, demo, pops = None):
  samples = sorted(list(tree.samples()))
  nodes_times = iter(sorted([(node, round(tree.time(node),5)) for node in tree.nodes() if tree.children(node)]))
  origin = Bundle(None, samples, pops)
  
  # backward in time
  bundle = origin
  node, time = next(nodes_times)
  for phase in demo.phases:
    bundle = bundle.transit(phase)
    while time < phase.t_end:
      bundle.coalesce(time, tree.children(node), node)
      try: node, time = next(nodes_times)
      except StopIteration: break
  
  # forward in time
  root = bundle
  root.root()
  while bundle.phase is not None:
    bundle.evolve()
    bundle.migrate()
    bundle = bundle.child
  
  # backward in time
  origin.origin()
  while bundle.parent:
    bundle = bundle.parent
    bundle.get_logp()
  
  return root.logp

def glike_trees(trees, demo, pops = None): # trees: generator or list of trees
  logp = 0
  for tree in trees:
    logp += glike(tree, demo, pops = pops)
  
  return logp


import msprime
demography = twoway_admixture_demography(20, 500, 0.3, 2000, 10000, 30000, 5000)
arg = msprime.sim_ancestry({"O": 1000, "A":1000, "B":1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
trees = [arg.at(pos).copy() for pos in range(int(3e5), int(3e7), int(3e5))]
pops = [0]*1000 + [1]*1000 + [2]*1000

demo = twoway_admixture_demo(20, 500, 0.3, 2000, 10000, 30000, 5000)
glike_trees(trees, twoway_admixture_demo(20, 500, 0.3, 2000, 10000, 30000, 5000), pops = pops)

A = B = C = D = E = F = 0
