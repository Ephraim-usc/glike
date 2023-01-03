from .glike import *

"""
def floor(number, digits):
  u = 0.1**digits
  return round(math.floor(number/u) * u, digits)

def ceil(number, digits):
  u = 0.1**digits
  return round(math.ceil(number/u) * u, digits)
"""


class Search():
  def __init__(self, names, values, limits = None, names_fixed = None):
    self.names = names
    self.values = dict(zip(names, values))
    if limits is None:
      limits = [(0, math.inf) for _ in names]
    self.limits = dict(zip(names, limits))
    if names_fixed is None:
      names_fixed = []
    self.names_fixed = names_fixed
    self.lrs = {name:0.1 for name in self.names}
  
  def get(self):
    return list(self.values.values())
  
  def set(self, values):
    self.values = dict(zip(self.names, values))
  
  def limit(self, name):
    limit = self.limits[name]
    low = limit[0]; high = limit[1]
    if isinstance(low, str):
      low = self.values[low]
    if isinstance(high, str):
      high = self.values[high]
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
    self.lrs[name] = max(0.1, self.lrs[name] * 0.5)
  
  def all_slow(self):
    for name in self.names:
      if self.lrs[name] > 0.1:
        return False
    return True


def estimate(trees, model, search, pops = None, tolerance = 0):
  x = search.get()
  logp_max = glike_trees(trees, model(*x), pops, tolerance)
  print(str(x) + " " + str(logp_max), flush = True)
  
  x_prev = x.copy()
  for _ in range(100):
    for name in [name for name in search.names if name not in search.names_fixed]:
      x_up = search.up(name)
      logp_up = glike_trees(trees, model(*x_up), pops, tolerance)
      x_down = search.down(name)
      logp_down = glike_trees(trees, model(*x_down), pops, tolerance)
      
      if (logp_up > max(logp_down, logp_max)):
        search.set(x_up)
        search.faster(name)
        logp_max = logp_up
      elif (logp_down > max(logp_up, logp_max)):
        search.set(x_down)
        search.faster(name)
        logp_max = logp_down
      else:
        search.slower(name)
    
    x = search.get()
    print(str(x) + " " + str(logp_max), flush = True)
    
    if x_prev == x and search.all_slow():
      break
    x_prev = x.copy()
  
  return search.get(), logp_max
