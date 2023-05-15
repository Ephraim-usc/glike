import math
import itertools
import tskit
import numpy as np
import scipy
import scipy.special
import scipy.linalg

import npe # numpy C extension

# test if a numpy array is identity matrix
def is_identity(x):
  return (x.shape[0] == x.shape[1]) and (x == np.eye(x.shape[0])).all()

# wrapper around scipy.special.logsumexp while testing shortcut situations
def logsumexp(x):
  if len(x) == 0:
    return -math.inf
  elif len(x) == 1:
    return x[0]
  else:
    return scipy.special.logsumexp(x)

# probability of coalescence at time a, given coalescence rate n and growth rate gr
def logp_coal(n, gr, a):
  if gr == 0:
    return math.log(n)
  else:
    return math.log(n) + gr * a

# probability of non-coalescence between a and b, given coalescence rate n and growth rate gr
def logp_noncoal(n, gr, a, b):
  if gr == 0:
    return - n * (b - a)
  else:
    return - n/gr * (math.exp(b*gr) - math.exp(a*gr))


class Phase:
  # t: beginning time of Phase
  # ns: list or 1d array of coalescent rates
  # grs: list or 1d array of growth rates
  # P: transition matrix at beginning of Phase
  # Q: continuous migration matrix, will be discretized when added into demo
  # populations: names of populations
  def __init__(self, t, ns, P = None, Q = None, populations = None):
    self.t = t
    self.t_end = math.inf
    self.parent = None
    self.child = None
    
    self.ns = np.array(ns)
    self.K = len(ns) # number of populations (during continuous period)
    
    if grs is not None:
      if len(grs) != self.K:
        raise Exception("Cannot initialize phase: len(grs) should equal len(ns)!!")
      self.grs = np.array(grs)
    else:
      self.grs = np.zeros(self.K)
    
    if P is not None:
      if P.shape[1] != self.K:
        raise Exception("Cannot initialize phase: P.shape[1] should equal len(ns)!")
      self.P = P
    else:
      self.P = np.identity(self.K)
    
    with np.errstate(divide='ignore'): # log(0) = -inf is expected
      self.logP = np.log(self.P)
    
    if Q is not None:
      if Q.shape[0] != self.K or Q.shape[1] != self.K
        raise Exception("Cannot initialize phase: Q.shape[0] and Q.shape[1] should both equal len(ns)!")
      self.Q = Q
    else:
      self.Q = None
    
    if populations is not None:
      assert len(populations) == self.K
      self.populations = populations
    else:
      self.populations = ['ABCDEFGHIJKLMNOPQRSTUVWXYZ'[i] for i in range(self.K)]
  
  # this function computes the logp of a genealogical tree during a non-migration period, given populations of lineages
  # N: number of lineages at t_end
  def logp(self, pop, N, coals):
    n = self.ns[pop]
    gr = self.grs[pop]
    coals.sort(reverse=True)
    
    buffer = 0
    b = self.t_end
    for a, incr in coals:
      buffer += N*(N-1)/2 * logp_noncoal(n, gr, a - self.t, b - self.t) if (b > a and N >= 2) else 0
      buffer += incr * logp_coal(n, gr, a) if incr >= 1 else 0
      N += incr; b = a
    buffer += N*(N-1)/2 * logp_noncoal(n, gr, 0, b - self.t) if N >= 2 else 0
    return buffer

class Demo:
  def __init__(self, demography = None):
    self.phases = list()
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
    ns = list()
    populations = list()
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
    self.logp = math.nan # evolution (coalescence and non-coalescence) log probability of this State
    self.logu = math.nan # absolute probability P(state|origin)
    self.logv = math.nan # absolute probability P(roots|state)
    self.children = list()
    self.parents = list()

class Bundle:
  def __init__(self, phase):
    self.phase = phase
    self.t = phase.t
    self.t_end = phase.t_end
    
    self.N = 0
    self.dict = dict() # for fast finding index of lineage
    self.lins = list() # list of lineages
    self.dests = list() # list of sets of descendents for each subtree
    self.ghosts = list() # indices of ghost lineages that are going to be deleted
    
    self.coals = list() # list of lists of (time, num_children - 1) tuples of nodes for each subtree
    self.mask = np.zeros([0, self.phase.K], dtype = np.int8) # N x K binary mask matrix of possible populations
    self.logmask = np.log(self.mask)
    
    self.states = dict() # dict of (value, state) tuples
    self.parent = None
    self.child = None
  
  def refresh(self):
    if self.ghosts:
      ghosts = set(self.ghosts); self.ghosts = list()
      self.lins = [lin for i, lin in enumerate(self.lins) if i not in ghosts]
      self.dests = [dests for i, dests in enumerate(self.dests) if i not in ghosts]
      self.coals = [coal for i, coal in enumerate(self.coals) if i not in ghosts]
      self.mask = np.array([mask_ for i,mask_ in enumerate(self.mask) if i not in ghosts])
    
    self.N = len(self.lins)
    self.dict = {lin:i for i, lin in enumerate(self.lins)}
    with np.errstate(divide='ignore'):
      self.logmask = np.log(self.mask)
  
  # backward in time through mass migrations
  # creating new Bundle, computing possible populations for each lineage
  def transit(self):
    bundle = Bundle(self.phase.parent)
    bundle.lins = self.lins.copy(); bundle.refresh()
    bundle.dests = [{lin} for lin in bundle.lins]
    bundle.coals = [list() for _ in bundle.lins]
    bundle.mask = (np.dot(self.mask, self.phase.parent.P) > 0).astype(np.int8)
    
    bundle.child = self
    self.parent = bundle
  
  # backward in time through continuous demography
  # summarizing coalescent results
  def coal(self, time, parent, children, pop):
    self.dict[parent] = len(self.lins)
    self.lins.append(parent)
    
    if children:
      idx = [self.dict[child] for child in children]
      self.ghosts.extend(idx)
      self.dests.append(set().union(*[self.dests[i] for i in idx]))
      self.coals.append([coal for i in idx for coal in self.coals[i]] + [(time, len(children) - 1)])
      mask_ = self.mask[idx,:].prod(axis = 0, dtype = np.int8)
    else:
      self.dests.append(set())
      self.coals.append([(time, - 1)])
      mask_ = np.ones(self.phase.K, dtype = np.int8)
    
    if pop != None: # designated by user
      if type(pop) is str:
        pop = self.phase.populations.index(pop)
      mask_[np.arange(len(mask_)) != pop] = 0 # set all other elements to zero, except the pop-th
    
    self.mask = np.vstack([self.mask, mask_]) # this is slow!!!!!
  
  # make it the root bundle
  def root(self):
    self.states = dict()
    pops = [np.nonzero(x)[0] for x in self.mask]
    for value in itertools.product(*pops):
      state = State()
      self.states[value] = state
  
  # forward in time through continuous demography
  # computing coalescence and non-coalescence probabilities
  def evolve(self):
    K = self.phase.K
    for value, state in self.states.items():
      Ns = [0 for _ in range(K)]
      coals = [list() for _ in range(K)]
      for i in range(self.N):
        Ns[value[i]] += 1
        coals[value[i]].extend(self.coals[i])
      state.logp = 0
      for pop in range(K):
        state.logp += self.phase.logp(pop, Ns[pop], coals[pop])
  
  def emigrate(self):
    child = self.child
    
    self.num_links = 0.0
    for value, state in self.states.items():
      outs = [None for _ in child.lins]
      for i,pop in enumerate(value):
        for dest in self.dests[i]:
          idx = child.dict[dest]
          outs[idx] = pop
      
      state.logP = self.phase.logP.T[outs, :] + child.logmask
      state.num_links = (state.logP > -math.inf).sum(axis = 1).prod(dtype = float)
      self.num_links += state.num_links
    
    #print(f"[{self.t}~{self.t_end}gen {self.phase.K} populations] [{len(self.child.lins)}-{len(self.lins)} lineages] [{len(self.states)} states] [{np.format_float_scientific(self.num_links, precision=6)} links]", flush = True)
  
  def immigrate(self, MAX_LINKS = 1e5):
    if hasattr(self.phase, 'mode') and self.phase.mode == "deterministic":
      self.immigrate_deterministic()
    elif hasattr(self.phase, 'mode') and self.phase.mode == "stochastic":
      self.immigrate_stochastic(MAX_LINKS)
    else:
      parent = self.parent
      if parent.num_links <= MAX_LINKS:
        self.immigrate_deterministic()
      else:
        self.immigrate_stochastic(MAX_LINKS)
  
  def immigrate_deterministic(self):
    N = self.N
    parent = self.parent
    for _, state_parent in parent.states.items():
      values, logps = npe.product_det(state_parent.logP)
      logps = logps.sum(axis = 1)
      
      for value, logp in zip(values, logps):
        value = tuple(value)
        if value in self.states:
          state = self.states[value]
        else:
          state = State()
          self.states[value] = state
        state_parent.children.append((logp, state))
        state.parents.append((logp, state_parent))
  
  def immigrate_stochastic(self, MAX_LINKS):
    N = self.N
    parent = self.parent
    
    for _, state_parent in parent.states.items():
      W = np.exp(state_parent.logP) # can try different definitions
      w = W.sum(axis=1, keepdims=True)
      state_parent.W = W/w
      state_parent.logw = np.log(w).sum()
    parent.logw = logsumexp([state_parent.logw for state_parent in parent.states.values()])
    
    for _, state_parent in parent.states.items():
      num = np.random.binomial(MAX_LINKS, math.exp(state_parent.logw - parent.logw))
      if num == 0:
        continue
      
      values, ws = npe.product_sto(state_parent.W, num)
      values, index, counts = np.unique(values, return_index=True, return_counts=True, axis = 0) # slow.
      logws = np.log(ws).sum(axis = 1)[index]
      logps = state_parent.logP[np.arange(N)[:,None], values.T].sum(axis = 0)
      
      for value, count, logw, logp in zip(values, counts, logws, logps):
        value = tuple(value)
        if value in self.states:
          state = self.states[value]
        else:
          state = State()
          self.states[value] = state
        logp_adj = logp + math.log(count) - math.log(MAX_LINKS) - (state_parent.logw - parent.logw + logw) # last term is the log prob. of sampling this state
        state_parent.children.append((logp_adj, state))
        state.parents.append((logp_adj, state_parent))
  
  def evaluate_logv(self):
    if self.parent:
      for state in self.states.values():
        state.logv = logsumexp([logp + parent.logv for logp, parent in state.parents]) + state.logp
    else:
      for state in self.states.values():
        state.logv = state.logp
    self.logv = logsumexp([state.logv for state in self.states.values()])
  
  def evaluate_logu(self):
    if self.child:
      for state in self.states.values():
        state.logu = logsumexp([logp + child.logu for logp, child in state.children]) + state.logp
    else:
      for state in self.states.values():
        state.logu = state.logp
    self.logu = logsumexp([state.logu for state in self.states.values()])


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
    bundle.coal(t, node, tree.children(node), samples.get(node, None))
  bundle.refresh()
  root = bundle
  
  # forward in time
  bundle = root
  bundle.root()
  bundle.evolve()
  bundle.evaluate_logv()
  while bundle.child:
    bundle.emigrate()
    bundle = bundle.child
    bundle.immigrate()
    bundle.evolve()
    bundle.evaluate_logv()
  
  # backward in time
  bundle = origin
  while bundle:
    bundle.evaluate_logu()
    bundle = bundle.parent
  
  return root.logu

def glike_trees(trees, demo, samples = None, prune = 0): # trees: generator or list of trees
  logps = [glike(tree, demo, samples = samples) for tree in trees]
  logps.sort()
  logp = sum(logps[math.ceil(prune * len(logps)):])
  return logp


def glike_verbose(tree, demo, samples = None):
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
    bundle.coal(t, node, tree.children(node), samples.get(node, None))
  bundle.refresh()
  root = bundle
  
  # forward in time
  bundle = root
  bundle.root()
  bundle.evolve()
  bundle.evaluate_logv()
  while bundle.child:
    bundle.emigrate()
    print(f"{bundle.t}~{bundle.t_end}gen, {bundle.phase.K} populations, {len(bundle.child.lins)}-{len(bundle.lins)} lineages, {len(bundle.states)} states, {np.format_float_scientific(bundle.num_links, precision=6)} links.", flush = True)
    bundle = bundle.child
    bundle.immigrate()
    bundle.evolve()
    bundle.evaluate_logv()
  print(f"{bundle.t}~{bundle.t_end}gen, {bundle.phase.K} populations, {len(bundle.lins)} lineages, {len(bundle.states)} states.", flush = True)
  
  # backward in time
  bundle = origin
  while bundle:
    bundle.evaluate_logu()
    bundle = bundle.parent
  
  return root
