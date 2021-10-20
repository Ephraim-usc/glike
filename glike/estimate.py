from .glike import *

def estimate(trees, labels, lmp_generator, initial_parameters, epochs):
  n_parameters = len(initial_parameters)
  parameters = initial_parameters.copy()
  
  for epoch in range(epochs):
    for i in range(n_parameters):
      parameters_down = parameters.copy(); parameters_down[i] = parameters_down[i] * 0.9
      parameters_up = parameters.copy(); parameters_up[i] = parameters_up[i] * 1.1
      lmp = lmp_generator(*parameters)
      lmp_down = lmp_generator(*parameters_down)
      lmp_up = lmp_generator(*parameters_up)
      logP = loglike_trees(trees, labels, lmp, 1000, stop = epoch*2+2).mean()
      logP_down = loglike_trees(trees, labels, lmp_down, 1000, stop = epoch*2+2).mean()
      logP_up = loglike_trees(trees, labels, lmp_up, 1000, stop = epoch*2+2).mean()
      
      if logP_up > logP:
        parameters = parameters_up.copy()
      if logP_down > logP:
        parameters = parameters_down.copy()
      print(parameters)
  
  return parameters
