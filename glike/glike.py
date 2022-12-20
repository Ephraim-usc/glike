import math
import numpy as np
import tskit
import itertools
import scipy
import scipy.special

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
    return - scipy.integrate.quad(n, a, b)
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
    self.ghosts = [] # indices of ghost lineages that are going to be deleted
    
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
    self.dict[parent] = len(self.lins)
    self.lins.append(parent)
    
    idx = [self.dict[child] for child in children]
    self.dests.append(set.union(*[self.dests[i] for i in idx]))
    self.pops.append(set.intersection(*[self.pops[i] for i in idx]))
    self.coals.append([coal for i in idx for coal in self.coals[i]] + [t] * (len(children) - 1))
    
    self.N -= len(children) - 1
    self.ghosts.extend(idx)
  
  def simplify(self):
    for i in sorted(self.ghosts, reverse=True):
      del self.lins[i]
      del self.dests[i]
      del self.pops[i]
      del self.coals[i]
    self.dict = {lin:i for i, lin in enumerate(self.lins)}
    self.ghosts = []
  
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
    for value in itertools.product(*self.pops):
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
      for ins in itertools.product(*pops_child):
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
      state.logp = scipy.special.logsumexp(logps) + state.logp_evolve
    self.logp = scipy.special.logsumexp([state.logp for state in self.states.values()])


def glike(tree, demo, pops = None):
  samples = sorted(list(tree.samples()))
  nodes_times = iter(sorted([(node, round(tree.time(node),5)) for node in tree.nodes() if tree.children(node)]))
  origin = Bundle(None, samples, pops)
  
  # backward in time
  phases = iter(demo.phases)
  bundle = origin.transit(next(phases))
  for node, time in nodes_times:
    while time > bundle.phase.t_end:
      bundle.simplify()
      bundle = bundle.transit(next(phases))
    bundle.coalesce(time, tree.children(node), node)
  bundle.simplify()
  
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
