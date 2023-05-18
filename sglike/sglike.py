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
  def __init__(self, t, t_end, ns, grs = None, P = None, Q = None, populations = None):
    self.t = t
    self.t_end = t_end
    self.parent = None
    self.child = None
    
    self.ns = np.array(ns)
    self.K = len(ns) # number of populations (during continuous period)
    
    if grs is not None:
      if len(grs) != self.K:
        raise Exception("Cannot initialize phase: len(grs) should equal len(ns)!")
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
      if not (self.K == Q.shape[0] == Q.shape[1]):
        raise Exception("Cannot initialize phase: Q.shape[0] and Q.shape[1] should both equal len(ns)!")
      if Q.any():
        self.Q = Q
      else:
        self.Q = None
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
  def __init__(self):
    self.phases = list()
  
  def add_phase(self, phase, discretize = 100):
    if phase.Q is not None:
      ns, grs, P, Q, populations = phase.ns, phase.grs, phase.P, phase.Q, phase.populations
      self.add_phase(Phase(phase.t, min(phase.t + discretize, phase.t_end), ns, grs, P = np.dot(P, scipy.linalg.expm(Q * min(discretize, phase.t_end-phase.t))), populations = populations))
      for t in np.arange(phase.t + discretize, phase.t_end, discretize):
        self.add_phase(Phase(t, min(t + discretize, phase.t_end), ns*np.exp(grs*(t-phase.t)), grs, P = scipy.linalg.expm(Q * min(discretize, phase.t_end-t)), populations = populations))
      return
    
    if len(self.phases):
      child = self.phases[-1]
      if child.t_end != phase.t:
        raise Exception("Cannot add phase: t should equal to t_end of the last existing phase!")
      if child.K != phase.P.shape[0]:
        raise Exception("Cannot add phase: phase.P.shape[0] is not equal to the number of populations in the last existing phase!")
      child.t_end = phase.t
      child.parent = phase
      phase.child = child
    else:
      if phase.t != 0:
        raise Exception("Cannot add phase: the first phase should have zero t!")
      if not is_identity(phase.P):
        raise Exception("Cannot add phase: the first phase should have identity P!")
    self.phases.append(phase)
  
  def print(self):
    import pandas as pd
    populations = None
    for phase in self.phases:
      print(f"[phase from {phase.t} to {phase.t_end}]")
      if not is_identity(phase.P):
        print(pd.DataFrame(phase.P, index = populations, columns = phase.populations))
      print(pd.DataFrame({"ns":phase.ns, "grs":phase.grs}, index = phase.populations).T)
      print("")
      populations = phase.populations


class State:
  def __init__(self):
    self.logp = math.nan # evolution (coalescence and non-coalescence) log probability of this State
    self.logu = math.nan # absolute probability P(state|origin), currently not used
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
    self.logmask = np.zeros([0, self.phase.K]) # N x K binary mask matrix of prior probabilities of population specification
    
    self.states = dict() # dict of (value, state) tuples
    self.parent = None
    self.child = None
  
  def refresh(self):
    if self.ghosts:
      ghosts = set(self.ghosts); self.ghosts = list()
      self.lins = [lin for i, lin in enumerate(self.lins) if i not in ghosts]
      self.dests = [dests for i, dests in enumerate(self.dests) if i not in ghosts]
      self.coals = [coal for i, coal in enumerate(self.coals) if i not in ghosts]
      self.logmask = np.array([logmask_ for i, logmask_ in enumerate(self.logmask) if i not in ghosts])
    
    self.N = len(self.lins)
    self.dict = {lin:i for i, lin in enumerate(self.lins)}
  
  # backward in time through mass migrations
  # creating new Bundle, computing possible populations for each lineage
  def transit(self):
    bundle = Bundle(self.phase.parent)
    bundle.lins = self.lins.copy(); bundle.refresh()
    bundle.dests = [{lin} for lin in bundle.lins]
    bundle.coals = [list() for _ in bundle.lins]
    
    with np.errstate(divide='ignore'):
      logmask_max = np.max(self.logmask, 1, keepdims=True)
      bundle.logmask = np.log(np.dot(np.exp(self.logmask - logmask_max), self.phase.parent.P)) + logmask_max
    
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
      logmask_ = self.logmask[idx,:].sum(axis = 0)
    else:
      self.dests.append(set())
      self.coals.append([(time, - 1)])
      logmask_ = np.zeros(self.phase.K)
    
    if pop != None: # designated by user
      if type(pop) is str:
        pop = self.phase.populations.index(pop)
      logmask_[np.arange(len(logmask_)) != pop] = -math.inf
    
    self.logmask = np.vstack([self.logmask, logmask_]) # this is slow
  
  # make it the root bundle
  def root(self):
    self.states = dict()
    pops = [np.where(x > -math.inf)[0] for x in self.logmask]
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
      state.logP = self.phase.logP.T[outs, :]
      self.num_links += (state.logP + child.logmask > -math.inf).sum(axis = 1).prod(dtype = float)
  
  def immigrate(self, MAX_LINKS = 1e4):
    if self.parent.num_links <= MAX_LINKS:
      self.immigrate_deterministic()
    else:
      self.immigrate_stochastic(MAX_LINKS)
  
  def immigrate_deterministic(self):
    N = self.N
    parent = self.parent
    for _, state_parent in parent.states.items():
      logPm = state_parent.logP.copy(); logPm[self.logmask == -math.inf] = -math.inf # masked logP
      values, logps = npe.product_det(logPm)
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
    
    minlogmask = math.log(1e-5 * self.t_end)
    logmask = self.logmask.copy()
    logmask[np.logical_and(logmask < minlogmask, logmask > -math.inf)] = minlogmask
    
    for _, state_parent in parent.states.items():
      W = np.exp(state_parent.logP + logmask)
      w = W.sum(axis=1, keepdims=True)
      state_parent.W = W/w
    
    logvs = np.array([state_parent.logv for state_parent in parent.states.values()])
    tmp = np.exp(logvs - parent.logv); tmp /= tmp.sum()
    nums = np.random.multinomial(MAX_LINKS, tmp)
    
    for state_parent, num in zip(parent.states.values(), nums):
      if num == 0:
        continue
      
      values, ws, logps = npe.product_sto(state_parent.W, state_parent.logP, num)
      values, index, counts = np.unique(values, return_index=True, return_counts=True, axis = 0) # slow.
      logws = np.log(ws).sum(axis = 1)[index]
      logps = logps.sum(axis = 1)[index]
      
      for value, count, logw, logp in zip(values, counts, logws, logps):
        value = tuple(value)
        if value in self.states:
          state = self.states[value]
        else:
          state = State()
          self.states[value] = state
        logp_adj = logp + math.log(count) - math.log(MAX_LINKS) - (state_parent.logv - parent.logv + logw) # last term is the log prob. of sampling this state
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
  
  def print(self):
    print(f"{self.t_end}~{self.t}gen", flush = True)
    print(f"populations: {self.phase.populations}", flush = True)
    print(f"# of lineages: {len(self.lins)}~{sum([len(dests) for dests in self.dests])}", flush = True)
    print(f"# of states: {len(self.states)}", flush = True)
    print(f"logv: {self.logv}", flush = True)
    logvs = sorted([state.logv for state in self.states.values()], reverse = True)[:5]
    shares = [f"{math.exp(logv - self.logv)*100:.2f}%" for logv in logvs]
    print(f"shares of top 5 states: {shares}", flush = True)

def glike(tree, demo, samples = None, verbose = False):
  if samples is None:
    samples = {}
  
  if type(tree) != tskit.trees.Tree:
    raise Exception("glike input type error: tree should be of type tskit.trees.Tree!")
  if type(demo) != Demo:
    raise Exception("glike input type error: demo should be of type Demo!")
  if type(samples) != dict:
    raise Exception("glike input type error: samples should be of type dict!")
  
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
  if verbose:
      bundle.print()
  
  while bundle.child:
    bundle.emigrate()
    if verbose:
      print(f"{np.format_float_scientific(bundle.num_links, precision=6)} links\n", flush = True)
    bundle = bundle.child
    bundle.immigrate()
    bundle.evolve()
    bundle.evaluate_logv()
    if verbose:
      bundle.print()
  
  return origin.logv


def glike_trees(trees, demo, samples = None, prune = 0): # trees: generator or list of trees
  if type(prune) not in (int, float):
    raise Exception("glike_trees input type error: prune should be int or float!")
  if not 0 <= prune <= 1:
    raise Exception("glike_trees input error: prune should be within [0, 1]!")
  
  logps = [glike(tree, demo, samples = samples) for tree in trees]
  logps.sort()
  logp = sum(logps[math.ceil(prune * len(logps)):])
  return logp
