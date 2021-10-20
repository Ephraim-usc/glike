from .glike import *

def get_delta(epoch):
  buffer = 0.01 + 0.19 * math.exp(-epoch / 20)
  buffer = round(buffer, 4)
  return buffer

def get_ntrees(epoch):
  buffer = min(2 + epoch, 20)
  return buffer

def estimate(trees, labels, lmp_generator, initial_parameters, epochs):
  n_parameters = len(initial_parameters)
  parameters = initial_parameters.copy()
  
  for epoch in range(epochs):
    delta = get_delta(epoch)
    n_trees = get_ntrees(epoch)
    
    for i in range(n_parameters):
      if (parameters[i] > 0.5) and (parameters[i] < 1): # it's a proportion parameter within [0.5, 1]
        parameters_down = parameters.copy(); parameters_down[i] = 1 - (1-parameters_down[i]) * (1 + delta)
        parameters_up = parameters.copy(); parameters_up[i] = 1 - (1-parameters_up[i]) * (1 - delta)
      else:
        parameters_down = parameters.copy(); parameters_down[i] = parameters_down[i] * (1 - delta)
        parameters_up = parameters.copy(); parameters_up[i] = parameters_up[i] * (1 + delta)
      
      lmp = lmp_generator(*parameters)
      lmp_down = lmp_generator(*parameters_down)
      lmp_up = lmp_generator(*parameters_up)
      logP = loglike_trees(trees, labels, lmp, 1000, stop = n_trees).mean()
      logP_down = loglike_trees(trees, labels, lmp_down, 1000, stop = n_trees).mean()
      logP_up = loglike_trees(trees, labels, lmp_up, 1000, stop = n_trees).mean()
      
      if logP_up > max(logP, logP_down):
        parameters = parameters_up.copy()
      if logP_down > max(logP, logP_up):
        parameters = parameters_down.copy()
      print(str(logP_down) + " " + str(logP) + " " + str(logP_up), flush = True)
      print(parameters, flush = True)
  
  return parameters
