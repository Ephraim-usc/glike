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

def generate_up_down(values, limits, precisions, pos):
  values_up = values.copy(); values_down = values.copy()
  values_up[pos] = min(limits[pos][1], values[pos] + precisions[pos])
  values_down[pos] = max(limits[pos][0], values[pos] - precisions[pos])
  return values_up, values_down

def estimate(trees, model, samples, transform, limits, precisions, flow = 10000, spread = 1e-5, prune = 0.5, verbose = False):
  _population_size = 10
  _num_offsprings = 50
  _num_epochs = 10
  population = []; scores = []
  
  print("Estimating parameters. Step 1: generating initial population (Genetic Algorithm)", flush = True)
  for _ in range(_population_size):
    values = [random.uniform(low, high) for low, high in limits]
    logp = glike_trees(trees, model(*transform(values)), samples = samples, flow = flow, spread = spread, prune = prune)
    if verbose:
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
    if verbose:
      print(f"{round_sig(transform(values))}, {logp}", flush = True)
    
    if logp > scores[-1]:
      population[-1] = values; scores[-1] = logp
      population = [genome for score, genome in sorted(zip(scores, population), reverse = True)]
      scores = [score for score, genome in sorted(zip(scores, population), reverse = True)]
    
    if _%_population_size == _population_size - 1:
      print(f"\nPopulation at generation {(_+1)//_population_size}:", flush = True)
      for genome, score in zip(population, scores):
        print(f"{round_sig(transform(genome))}, {score}", flush = True)
      print("\n")
  
  print("Estimating parameters. Step 3: Refining best genome", flush = True)
  values = population[0]; score = scores[0]
  for _ in range(_num_epochs):
    for pos in range(len(values)):
      values_up, values_down = generate_up_down(values, limits, precisions, pos)
      logp_up = glike_trees(trees, model(*transform(values_up)), samples = samples, flow = flow, spread = spread, prune = prune)
      if verbose:
        print(f"{round_sig(transform(values_up))}, {logp_up}", flush = True)
      logp = glike_trees(trees, model(*transform(values)), samples = samples, flow = flow, spread = spread, prune = prune)
      if verbose:
        print(f"{round_sig(transform(values))}, {logp}", flush = True)
      logp_down = glike_trees(trees, model(*transform(values_down)), samples = samples, flow = flow, spread = spread, prune = prune)
      if verbose:
        print(f"{round_sig(transform(values_down))}, {logp_down}", flush = True)
      
      if (logp_up > max(logp_down, logp)):
        values = values_up.copy(); logp = logp_up
      elif (logp_down > max(logp_up, logp)):
        values = values_down.copy(); logp = logp_down
    
    print(f"\nCurrent best at epoch {_+1}:", flush = True)
    print(f"{round_sig(transform(values))}, {logp}", flush = True)
    print("\n")
  
  return values, logp
