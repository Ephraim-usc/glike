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

# we do not use np.arange because it has errors at ~10 digits after the decimal point
def intervals(start, stop, step):
  a = []
  t = start
  while t < stop:
    a.append(t)
    t += step
  b = a[1:]+[stop]
  return list(zip(a,b))

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
    self.parent = None
    self.child = None
    
    if type(t) not in (int, float):
      raise Exception("Cannot initialize phase: t should be a number!")
    if type(t_end) not in (int, float):
      raise Exception("Cannot initialize phase: t_end should be a number!")
    if not (0 <= t < t_end):
      raise Exception("Cannot initialize phase: it is required that 0 <= t < t_end!")
    self.t = float(t)
    self.t_end = float(t_end)
    
    if type(ns) not in (list, np.ndarray):
      raise Exception("Cannot initialize phase: ns should be a list or a numpy array!")
    if (type(ns) is np.ndarray) and (ns.ndim != 1):
      raise Exception("Cannot initialize phase: ns should be a list or a 1d array!")
    if min(ns) <= 0:
      raise Exception("Cannot initialize phase: ns should be all positive!")
    self.ns = np.array(ns)
    self.K = len(ns) # number of populations (during continuous period)
    
    if grs is not None:
      if type(grs) not in (list, np.ndarray):
        raise Exception("Cannot initialize phase: grs should be a list or a numpy array!")
      if (type(grs) is np.ndarray) and (ns.ndim != 1):
        raise Exception("Cannot initialize phase: grs should be a list or a 1d array!")
      if min(grs) < 0:
        raise Exception("Cannot initialize phase: grs should be all non-negative!")
      if len(grs) != self.K:
        raise Exception("Cannot initialize phase: grs should be of equal length as ns!")
      self.grs = np.array(grs)
    else:
      self.grs = np.zeros(self.K)
    
    if P is not None:
      if (type(P) != np.ndarray) or (P.ndim != 2):
        raise Exception("Cannot initialize phase: P should be a 2d numpy array!")
      if P.shape[1] != self.K:
        raise Exception("Cannot initialize phase: P.shape[1] should equal len(ns)!")
      if P.min() < -1e-8: # as P is in numpy array format, we allow small numerical errors
        raise Exception("Cannot initialize phase: P should be all non-negative!")
      self.P = P
    else:
      self.P = np.identity(self.K)
    
    with np.errstate(divide='ignore'): # log(0) = -inf is expected
      self.logP = np.log(self.P)
    
    if Q is not None:
      if (type(Q) != np.ndarray) or (Q.ndim != 2):
        raise Exception("Cannot initialize phase: Q should be a 2d numpy array!")
      if not (self.K == Q.shape[0] == Q.shape[1]):
        raise Exception("Cannot initialize phase: Q.shape[0] and Q.shape[1] should both equal len(ns)!")
      if Q.any():
        self.Q = Q
      else:
        self.Q = None
    else:
      self.Q = None
    
    if populations is not None:
      if type(populations) not in (list, np.ndarray):
        raise Exception("Cannot initialize phase: populations should be a list or a numpy array!")
      if (type(populations) is np.ndarray) and (populations.ndim != 1):
        raise Exception("Cannot initialize phase: populations should be a list or a 1d array!")
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
      ts = intervals(phase.t, phase.t_end, discretize); t, t_end = ts[0]
      self.add_phase(Phase(t, t_end, ns, grs, P = np.dot(P, scipy.linalg.expm(Q * (t_end - t))), populations = populations))
      for t, t_end in ts[1:]:
        self.add_phase(Phase(t, t_end, ns*np.exp(grs*(t-phase.t)), grs, P = scipy.linalg.expm(Q * (t_end - t)), populations = populations))
      return
    
    if len(self.phases):
      child = self.phases[-1]
      if child.t_end != phase.t:
        raise Exception("Cannot add phase: t should equal to t_end of the last existing phase!")
      if child.K != phase.P.shape[0]:
        raise Exception("Cannot add phase: phase.P.shape[0] should equal the number of populations of the last existing phase!")
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
      logmask_max[logmask_max == -np.inf] = 0 # this is to avoid the (-inf - -inf) error
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
  
  def immigrate(self, flow, spread):
    if self.parent.num_links <= flow:
      self.immigrate_deterministic()
    else:
      self.immigrate_stochastic(flow, spread)
  
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
  
  def immigrate_stochastic(self, flow, spread):
    N = self.N
    parent = self.parent
    
    minlogmask = math.log(spread * self.t_end)
    logmask = self.logmask.copy()
    logmask[np.logical_and(logmask < minlogmask, logmask > -math.inf)] = minlogmask
    
    for _, state_parent in parent.states.items():
      W = np.exp(state_parent.logP + logmask)
      w = W.sum(axis=1, keepdims=True)
      state_parent.W = W/w
    
    logvs = np.array([state_parent.logv for state_parent in parent.states.values()])
    tmp = np.exp(logvs - parent.logv); tmp /= tmp.sum()
    nums = np.random.multinomial(flow, tmp)
    
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
        logp_adj = logp + math.log(count) - math.log(flow) - (state_parent.logv - parent.logv + logw) # last term is the log prob. of sampling this state
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


def glike(tree, demo, samples = None, flow = 1e4, spread = 1e-5, verbose = False):
  if samples is None:
    samples = {}
  
  if type(tree) != tskit.trees.Tree:
    raise Exception("glike input type error: tree should be of type tskit.trees.Tree!")
  if type(demo) != Demo:
    raise Exception("glike input type error: demo should be of type Demo!")
  if type(samples) != dict:
    raise Exception("glike input type error: samples should be of type dict!")
  if type(flow) not in (int, float):
    raise Exception("glike input type error: flow should be an int!")
  if (type(spread) not in (int, float)) or (spread > 1):
    raise Exception("glike input type error: spread should be a number between 0 and 1 (e.g., 1e-5)!")
  
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
    bundle.immigrate(flow, spread)
    bundle.evolve()
    bundle.evaluate_logv()
    if verbose:
      bundle.print()
  
  return origin.logv

def glike_trees(trees, demo, samples = None, flow = 1e4, spread = 1e-5, prune = 0): # trees: generator or list of trees
  if type(prune) not in (int, float):
    raise Exception("glike_trees input type error: prune should be int or float!")
  if not 0 <= prune <= 1:
    raise Exception("glike_trees input error: prune should be within [0, 1]!")
  
  logps = [glike(tree, demo, samples = samples, flow = flow, spread = spread) for tree in trees]
  logps.sort()
  logp = sum(logps[math.ceil(prune * len(logps)):])
  return logp
