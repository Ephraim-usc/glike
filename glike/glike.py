import math
import numpy as np
import tskit
import itertools
import scipy
import scipy.special

def is_identity(x):
  return (x.shape[0] == x.shape[1]) and (x == np.eye(x.shape[0])).all()

# probability of coalescence at time a
def logp_n(n, a):
  if type(n) is float:
    return math.log(n)
  elif type(n) is tuple:
    return math.log(n[0]) + n[1] * a
  elif callable(n):
    return math.log(n(a))
  else:
    raise Exception("not supported n type: not a float, a tuple or a callable")

# probability of non-coalescence between a and b
def logp_intn(n, a, b):
  if type(n) is float:
    return - n * (b - a)
  if type(n) is tuple:
    return - n[0] / n[1] * (math.exp(b*n[1]) - math.exp(a*n[1]))
  elif callable(n):
    return - scipy.integrate.quad(n, a, b)
  else:
    raise Exception("not supported n type: not a float, a tuple or a callable")


class Phase:
  # t: beginning time of Phase
  # ns: list of coalescent rates, each could be a number, a (initial_rate, growth_rate) tuple, or any function f(t)
  # P: transition matrix at beginning of Phase
  def __init__(self, t, ns, P = None, populations = None):
    self.t = t
    self.t_end = math.inf
    self.parent = None
    self.child = None
    
    self.ns = ns
    self.K = len(ns) # number of populations (during continuous period)
    
    if P is not None:
      assert P.shape[1] == self.K
      self.P = P
    else:
      self.P = np.identity(self.K)
    
    with np.errstate(divide='ignore'):
      self.ins_outs_mapping()
      logP = np.log(self.P)
      logP[np.isneginf(logP)] = 0 # this is necessary because otherwise "-inf * 0 = nan" will cause trouble
      self.logP = logP
    
    if populations is not None:
      assert len(populations) == self.K
      self.populations = populations
    else:
      self.populations = ['ABCDEFGHIJKLMNOPQRSTUVWXYZ'[i] for i in range(self.K)]
  
  # build fasting mappings of non-zero migrations in P
  def ins_outs_mapping(self):
    self.ins2outs = [set() for _ in range(self.P.shape[0])]
    self.outs2ins = [set() for _ in range(self.P.shape[1])]
    ins, outs = np.where(self.P > 0)
    for in_, out in zip(ins, outs):
      self.ins2outs[in_].add(out)
      self.outs2ins[out].add(in_)
  
  # this function computes the logp of a genealogical tree during a non-migration period, given populations of lineages
  # N: number of lineages at t_end
  def logq(self, pop, N, births):
    n = self.ns[pop]
    births.sort(reverse=True)
    
    buffer = 0
    time_ = self.t_end
    for time, delta in births:
      buffer += 0.5*N*(N-1) * logp_intn(n, time, time_) if time_ > time and N >= 2 else 0 # think of a better way to deal with present samples!!!
      buffer += delta * logp_n(n, time) if delta >= 1 else 0
      N += delta; time_ = time
    buffer += 0.5*N*(N-1) * logp_intn(n, self.t, time_) if N >= 2 else 0
    return buffer

class Demo:
  def __init__(self, demography = None):
    self.phases = []
    if demography:
      self.from_demography(demography)
  
  def add_phase(self, phase):
    if len(self.phases):
      child = self.phases[-1]
      assert child.t < phase.t, "time error when adding phase!"
      assert child.K == phase.P.shape[0], "shape error when adding phase!"
      child.t_end = phase.t
      child.parent = phase
      phase.child = child
    else:
      assert is_identity(phase.P), "first phase should have identity P!"
    self.phases.append(phase)
  
  def at(self, t):
    for phase in self.phases:
      if t < phase.t_end:
        return phase
  
  def print(self):
    import pandas as pd
    populations = None
    for phase in self.phases:
      print(f"phase at {phase.t}")
      if not is_identity(phase.P):
        print(pd.DataFrame(phase.P, index = populations, columns = phase.populations))
      print(dict(zip(phase.populations, phase.ns)))
      print("")
      populations = phase.populations
  
  def from_demography(self, demography):
    import msprime
    import pandas as pd
    t = 0
    ns = []
    populations = []
    for population in demography.populations:
      populations.append(population.name)
      if population.growth_rate != 0:
        ns.append((1/population.initial_size, population.growth_rate))
      else:
        ns.append(1/population.initial_size)
    self.add_phase(Phase(0, ns, populations = populations))
    
    for event in demography.events:
      if type(event) is msprime.demography.Admixture:
        populations_ = populations; populations = populations.copy(); 
        populations.remove(event.derived)
        P = pd.DataFrame(0, index = populations_, columns = populations)
        for population in populations:
          P.loc[population, population] = 1
        for anc, prop in zip(event.ancestral, event.proportions):
          P.loc[event.derived, anc] = prop
        P = P.values
      elif type(event) is msprime.demography.PopulationSplit:
        populations_ = populations; populations = populations.copy()
        for derived in event.derived:
          populations.remove(derived)
        P = pd.DataFrame(0, index = populations_, columns = populations)
        for population in populations:
          P.loc[population, population] = 1
        for derived in event.derived:
          P.loc[derived, event.ancestral] = 1
        P = P.values
      elif type(event) is msprime.demography.PopulationParametersChange:
        populations_ = populations; populations = populations.copy()
        P = None
      else:
        continue
      
      t_ = t; t = event.time
      ns_ = ns.copy()
      ns = [None for _ in populations]
      for population in populations:
        n = ns_[populations_.index(population)]
        if type(n) is float:
          ns[populations.index(population)] = n
        else:
          ns[populations.index(population)] = (math.log(n[0]) + n[1] * (t-t_), n[1])
      
      if type(event) is msprime.demography.PopulationParametersChange:
        if event.population == -1:
          tmp = populations.copy()
        else:
          tmp = [event.population]
        for population in tmp:
          ns[populations.index(population)] = (1/event.initial_size, event.growth_rate) if event.growth_rate and event.growth_rate>0 else 1/event.initial_size
      
      self.add_phase(Phase(t, ns, P = P, populations = populations))


class State:
  def __init__(self):
    self.logq = math.nan # evolution (coalescence and non-coalescence) log probability of this State
    self.logp = math.nan # absolute probability P(state|origin)
    self.logp_ = math.nan # absolute probability P(root|state), currently not used
    self.children = []

class Bundle:
  # phase: the Phase this Bundle is in, the origin Bundle should have None phase
  # lins: list of lineages
  # pops: list of deterministic populations for each lineage
  def __init__(self, phase):
    self.phase = phase
    self.t = phase.t
    self.t_end = phase.t_end
    
    self.N = 0
    self.dict = dict() # for fast finding index of lineage
    self.lins = [] # list of lineages
    self.pops = [] # list of sets of possible populations for each subtree
    self.dests = [] # list of sets of descendents for each subtree
    self.births = [] # list of lists of (time, num_children - 1) tuples of nodes for each subtree
    
    self.ghosts = [] # indices of ghost lineages that are going to be deleted
    self.states = {} # dict of (value, state) tuples
    self.parent = None
    self.child = None
  
  def refresh(self):
    if self.ghosts:
      ghosts = set(self.ghosts); self.ghosts = []
      self.lins = [lin for i, lin in enumerate(self.lins) if i not in ghosts]
      self.dests = [dests for i, dests in enumerate(self.dests) if i not in ghosts]
      self.pops = [pops for i, pops in enumerate(self.pops) if i not in ghosts]
      self.births = [births for i, births in enumerate(self.births) if i not in ghosts]
    
    self.N = len(self.lins)
    self.dict = {lin:i for i, lin in enumerate(self.lins)}
  
  # backward in time through mass migrations
  # creating new Bundle, computing possible populations for each lineage
  def transit(self):
    bundle = Bundle(self.phase.parent)
    bundle.lins = self.lins.copy(); bundle.refresh()
    bundle.dests = [{lin} for lin in bundle.lins]
    bundle.births = [[] for _ in bundle.lins]
    bundle.pops = [set().union(*[bundle.phase.ins2outs[in_] for in_ in pops]) for pops in self.pops] # set().union() can handle empty pops
    
    bundle.child = self
    self.parent = bundle
  
  # backward in time through continuous demography
  # summarizing coalescent results
  def birth(self, time, parent, children, pop):
    self.dict[parent] = len(self.lins)
    self.lins.append(parent)
    
    idx = [self.dict[child] for child in children]
    self.ghosts.extend(idx)
    self.dests.append(set().union(*[self.dests[i] for i in idx]))
    self.births.append([birth for i in idx for birth in self.births[i]] + [(time, len(children) - 1)])
    
    if pop: # designated by user
      if type(pop) is str:
        pop = self.phase.populations.index(pop)
      self.pops.append({pop})
    elif children: # intersection of children pops
      self.pops.append(set.intersection(*[self.pops[i] for i in idx]))
    else: # end node without designation, could be in any population
      self.pops.append(set(range(len(self.phase.ns))))
  
  # make it the root bundle
  def root(self):
    for value in itertools.product(*self.pops):
      state = State()
      self.states[value] = state
  
  # forward in time through continuous demography
  # computing coalescence and non-coalescence probabilities
  def evolve(self):
    K = self.phase.K
    for value, state in self.states.items():
      Ns = [0 for _ in range(K)]
      births = [[] for _ in range(K)]
      for i in range(self.N):
        Ns[value[i]] += 1
        births[value[i]].extend(self.births[i])
      state.logq = 0
      for pop in range(K):
        state.logq += self.phase.logq(pop, Ns[pop], births[pop])
  
  # forward in time through mass migrations
  # creating children states and computing migration probabilities
  def migrate(self):
    outs2ins = self.phase.outs2ins
    logP = self.phase.logP
    child = self.child
    for value, state in self.states.items():
      outs = [None for _ in child.lins]
      for i,pop in enumerate(value):
        for dest in self.dests[i]:
          outs[child.dict[dest]] = pop
      
      pops_child = [set.intersection(child.pops[i], outs2ins[outs[i]]) for i in range(child.N)]
      for ins in itertools.product(*pops_child):
        migrations = np.zeros(self.phase.P.shape)
        for in_, out in zip(ins, outs):
          migrations[in_, out] += 1
        logq = (logP * migrations).sum()
        
        state_child = child.states.setdefault(ins, State())
        state.children.append((logq, state_child))
  
  def evaluate(self):
    if not self.child:
      for state in self.states.values():
        state.logp = state.logq
    else:
      for state in self.states.values():
        state.logp = scipy.special.logsumexp([logq + child.logp for logq, child in state.children]) + state.logq
    
    if not self.parent:
      self.logp = scipy.special.logsumexp([state.logp for state in self.states.values()]) if self.states else -math.inf


def glike(tree, demo, samples = None):
  if not samples:
    samples = {}
  times_nodes = iter(sorted([(round(tree.time(node),5), node) for node in tree.nodes()]))
  origin = Bundle(demo.phases[0])
  
  # backward in time
  bundle = origin
  for t, node in times_nodes:
    while t >= bundle.phase.t_end:
      bundle.refresh()
      bundle.transit()
      bundle = bundle.parent
    bundle.birth(t, node, tree.children(node), samples.get(node, None))
  bundle.refresh()
  
  root = bundle
  root.root()
  if len(root.states) == 0:
    return -math.inf
  
  # forward in time
  bundle = root
  while bundle.child:
    bundle.migrate()
    bundle = bundle.child
    #print("bundle from {} to {} has {} states".format(bundle.t, bundle.t_end, len(bundle.states)), flush = True)
  
  # backward in time
  bundle = origin
  while bundle:
    bundle.evolve()
    bundle.evaluate()
    bundle = bundle.parent
  
  return root.logp

def glike_trees(trees, demo, samples = None, prune = 0): # trees: generator or list of trees
  logps = [glike(tree, demo, samples = samples) for tree in trees]
  logps.sort()
  logp = sum(logps[math.ceil(prune * len(logps)):])
  return logp

