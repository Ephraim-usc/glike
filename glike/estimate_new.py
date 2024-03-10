from .glike import *
import random

def round_sig(x, sig = 4):
  if type(x) is list:
    return [round_sig(x_, sig) for x_ in x]
  else:
    return round(x, sig-int(math.floor(math.log10(abs(x))))-1)

def generate_offspring(values1, values2, limits, precisions):
  values = []
  for value1, value2, limit, precision in zip(values1, values2, limits, precisions):
    value = value2 if np.random.rand() > 0.5 else value1
    value += precision * np.random.normal()
    value = max(value, limit[0]); value = min(value, limit[1])
    values.append(value)
  return values

def estimate(trees, model, samples, transform, limits, precisions, flow = 10000, spread = 1e-5, prune = 0.5, epochs = 100, verbose = False):
  _population_size = 10
  _num_offsprings = 100
  population = []; scores = []
  
  print("Estimating parameters. Step 1: generating initial population (Genetic Algorithm)", flush = True)
  for _ in range(_population_size):
    values = [random.uniform(low, high) for low, high in limits]
    logp = glike_trees(trees, model(*transform(values)), samples = samples, flow = flow, spread = spread, prune = prune)
    print(f"{round_sig(transform(values))}, {logp}", flush = True)
    population.append(values); scores.append(logp)
  
  population = [genome for score, genome in sorted(zip(scores, population), reverse = True)]
  scores = [score for score, genome in sorted(zip(scores, population), reverse = True)]
  
  print("\nInitial population:", flush = True)
  for genome, score in zip(population, scores):
    print(f"{round_sig(transform(genome))}, {score}", flush = True)
  print("\n")
  
  print("Estimating parameters. Step 2: generating offsprings (Genetic Algorithm)", flush = True)
  for _ in range(_num_offsprings):
    idx1, idx2 = random.randrange(_population_size), random.randrange(_population_size)
    values1, values2 = population[idx1], population[idx2]
    values = generate_offspring(values1, values2, limits, precisions)
    logp = glike_trees(trees, model(*transform(values)), samples = samples, flow = flow, spread = spread, prune = prune)
    print(f"{round_sig(transform(values))}, {logp}", flush = True)
    
    if logp > scores[-1]:
      population[-1] = values; scores[-1] = logp
      population = [genome for score, genome in sorted(zip(scores, population), reverse = True)]
      scores = [score for score, genome in sorted(zip(scores, population), reverse = True)]
    
    if _ > 0 and _%_population_size == 0:
      print("\nPopulation at generation {_//_population_size}:", flush = True)
      for genome, score in zip(population, scores):
        print(f"{round_sig(transform(genome))}, {score}", flush = True)
      print("\n")
  
  print("Estimating parameters. Step 3: Refining best genome", flush = True)
  



class Search():
  def __init__(self, names, values, limits = None, names_fixed = None, precision = 0.05):
    self.names = names
    self.values = dict(zip(names, values))
    if limits is None:
      limits = [(0, math.inf) for _ in names]
    self.limits = dict(zip(names, limits))
    if names_fixed is None:
      names_fixed = []
    self.names_fixed = names_fixed
    self.lrs = {name:0.1 for name in self.names}
    self.precision = precision
  
  def get(self):
    return list(self.values.values())
  
  def set(self, values):
    self.values = dict(zip(self.names, values))
  
  def limit(self, name):
    limit = self.limits[name]
    low = limit[0]; high = limit[1]
    if isinstance(low, str):
      low = eval(low, self.values.copy())
    if isinstance(high, str):
      high = eval(high, self.values.copy())
    return low, high
  
  def up(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = round(min(high - 1e-5, value + step), 5)
    return list(values.values())
  
  def down(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      step = (value - low) * lr
    else:
      step = (high - value) * lr
    step = max(step, 1e-5)
    values = self.values.copy()
    values[name] = round(max(low + 1e-5, value - step), 5)
    return list(values.values())
  
  def faster(self, name):
    self.lrs[name] = min(0.5, self.lrs[name] * 1.5)
  
  def slower(self, name):
    self.lrs[name] = max(self.precision, self.lrs[name] * 0.5)
  
  def cold(self):
    for name in self.names:
      if (name not in self.names_fixed) and (self.lrs[name] > self.precision):
        return False
    return True


def estimate(trees, model, search, samples = None, flow = 10000, spread = 1e-5, prune = 0.5, epochs = 100, verbose = False):
  x = search.get()
  logp = glike_trees(trees, model(*x), samples = samples, flow = flow, spread = spread, prune = prune)
  print(str(x) + " " + str(logp), flush = True)
  
  xs = []
  logps = []
  for _ in range(epochs):
    for name in [name for name in search.names if name not in search.names_fixed]:
      x = search.get()
      logp = glike_trees(trees, model(*x), samples = samples, flow = flow, spread = spread, prune = prune)
      x_up = search.up(name)
      logp_up = glike_trees(trees, model(*x_up), samples = samples, flow = flow, spread = spread, prune = prune)
      x_down = search.down(name)
      logp_down = glike_trees(trees, model(*x_down), samples = samples, flow = flow, spread = spread, prune = prune)
      
      if verbose:
        print(" ", flush = True)
        print("x_up: " + str(x_up) + " " + str(logp_up), flush = True)
        print("x: " + str(x) + " " + str(logp), flush = True)
        print("x_down: " + str(x_down) + " " + str(logp_down), flush = True)
        print(" ", flush = True)
      
      if (logp_up > max(logp_down, logp)):
        search.set(x_up)
        search.faster(name)
      elif (logp_down > max(logp_up, logp)):
        search.set(x_down)
        search.faster(name)
      else:
        search.slower(name)
    
    x = search.get()
    logp = glike_trees(trees, model(*x), samples = samples, flow = flow, spread = spread, prune = prune)
    xs.append(x); logps.append(logp)
    print(str(x) + " " + str(logp), flush = True)
    
    if len(logps) >= 5 and sum(logps[-5:-3]) >= sum(logps[-2:]):
      break
  
  idx = logps.index(max(logps))
  x, logp = xs[idx], logps[idx]
  return x, logp
