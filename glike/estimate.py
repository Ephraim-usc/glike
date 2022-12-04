from .glike import *


class Searchspace():
  def __init__(self, names, values, limits = None):
    self.names = names
    self.values = dict(zip(names, values))
    if limits is None:
      limits = [(0, math.inf) for _ in names]
    self.limits = dict(zip(names, limits))
    self.lrs = {name:0.1 for name in self.names}
  
  def now(self):
    return list(self.values.values())
  
  def assign(self, values):
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
      value = low + (value - low) * (1 + lr)
    else:
      value = high - (high - value) * (1 - lr)
    values = self.values.copy()
    values[name] = round(value, 6)
    return list(values.values())
  
  def down(self, name):
    value = self.values[name]
    lr = self.lrs[name]
    low, high = self.limit(name)
    if value < (low + high)/2:
      value = low + (value - low) * (1 - lr)
    else:
      value = high - (high - value) * (1 + lr)
    values = self.values.copy()
    values[name] = round(value, 6)
    return list(values.values())
  
  def faster(self, name):
    self.lrs[name] = min(0.5, self.lrs[name] * 1.5)
  
  def slower(self, name):
    self.lrs[name] = max(0.1, self.lrs[name] * 0.5)



def estimate(model, trees, x_ref, x, fixed):
  steps = [new_step_func("regular") for _ in range(2)] + [new_step_func("porportion") for _ in range(2)] + [new_step_func("regular") for _ in range(len(x) - 4)]
  
  logp_ref = glike_trees(trees, model(*x_ref))
  print(str(x_ref) + " " + str(logp_ref))
  logp_max = glike_trees(trees, model(*x))
  print(str(x) + " " + str(logp_max))
  
  x_prev = x.copy()
  for _ in range(100):
    for i in [_ for _ in range(len(x)) if _ not in fixed]:
      x_inc = x.copy(); x_inc[i] = steps[i](x[i], True)
      logp_inc = glike_trees(trees, model(*x_inc))
      x_dec = x.copy(); x_dec[i] = steps[i](x[i], False)
      logp_dec = glike_trees(trees, model(*x_dec))
      
      if ((logp_inc > logp_dec) and (logp_inc > logp_max)):
        x = x_inc
        logp_max = logp_inc
        steps[i].lr = min(0.5, steps[i].lr * 1.5)
      if ((logp_dec > logp_inc) and (logp_dec > logp_max)):
        x = x_dec
        logp_max = logp_dec
        steps[i].lr = min(0.5, steps[i].lr * 1.5)
      if ((logp_max > logp_inc) and (logp_max > logp_dec)):
        pass
        steps[i].lr = max(0.1, steps[i].lr * 0.5)
    
    print(str(x) + " " + str(logp_max))
    #print([step.lr for step in steps])
    if x_prev == x and min([step.lr <= 0.1 for step in steps]):
      break
    x_prev = x

