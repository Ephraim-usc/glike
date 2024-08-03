from .glike import *

class Search():
  def __init__(self, x0, limits = None, names_fixed = None, precision = 0.05):
    self.names = list(x0.keys())
    self.values = x0.copy()
    if limits is None:
      limits = [(0, math.inf) for _ in self.names]
    self.limits = dict(zip(self.names, limits))
    if names_fixed is None:
      names_fixed = []
    self.names_fixed = names_fixed
    self.lrs = {name:0.1 for name in self.names}
    self.precision = precision
  
  def get(self):
    return self.values
  
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
    values[name] = min(high, round(value + step, 5))
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


def maximize(fun, x0, limits = None, precision = 0.05, epochs = 20, verbose = False):
  # fun: The objective function to be maximized.
  # x0: the dict of initial parameters, such that the initial output would be fun(**x0)
  # limits: the list of 2-tuples that defines the boundaries
  # precision: a float that defines the (proportional) step size
  # epochs: an integer that defines the maximum number of epochs
  # verbose: True if intermediate results are printed, False if not
  
  search = Search(x0, limits = limits, precision = precision)
  names = list(x0.keys())
  
  y0 = fun(**x0)
  print(str(x0) + " " + str(y0), flush = True)
  
  xs = []
  ys = []
  for _ in range(epochs):
    for name in names:
      x = search.get()
      y = fun(**x)
      x_up = search.up(name)
      y_up = fun(**x_up)
      x_down = search.down(name)
      y_down = fun(**x_down)
      
      if verbose:
        print(" ", flush = True)
        print("x_up: " + str(x_up) + " " + str(y_up), flush = True)
        print("x: " + str(x) + " " + str(y), flush = True)
        print("x_down: " + str(x_down) + " " + str(y_down), flush = True)
        print(" ", flush = True)
      
      if (y_up > max(y_down, y)):
        search.set(x_up)
        search.faster(name)
      elif (y_down > max(y_up, y)):
        search.set(x_down)
        search.faster(name)
      else:
        search.slower(name)
    
    x = search.get()
    y = fun(**x)
    xs.append(x); ys.append(y)
    print(str(x) + " " + str(y), flush = True)
    
    if len(ys) >= 5 and sum(ys[-5:-3]) >= sum(ys[-2:]):
      break
  
  idx = ys.index(max(ys))
  x, y = xs[idx], ys[idx]
  return x, y
