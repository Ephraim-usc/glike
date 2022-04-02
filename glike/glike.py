import math
import time
import numpy as np
import pandas as pd
import scipy.linalg
from scipy.special import logsumexp
import msprime
np.seterr(divide='ignore')




class Stage: #lineage Markov Process
  def __init__(self, in_pops, out_pops, P = None):
    self.in_pops = in_pops
    self.out_pops = out_pops
    
    if P is None:
      P = np.identity(len(in_pops))
    
    self.P = pd.DataFrame(P)
    self.P.index = in_pops; self.P.columns = out_pops
    
    self.Q = pd.DataFrame(np.zeros([len(out_pops), len(out_pops)]))
    self.Q.index = self.Q.columns = out_pops
    
    self.QQ = pd.DataFrame(np.zeros([len(out_pops)**2, len(out_pops)**2]))
    self.QQ.index = self.QQ.columns = [x+y for x in out_pops for y in out_pops]
  
  def add_coalescent_rate(self, pops, n):
    if type(pops) != list:
      pops = [pops]
    for pop in pops:
      for pop2 in pops:
      self.QQ.loc[pop+pop2, pop+pop2] = -n
  
  def add_migration_rate(self, source, dest, rate):
    self.Q.loc[source, dest] += rate
    self.Q.loc[source, source] -= rate
    
    for pop in self.out_pops:
      self.QQ.loc[]
  
  def add_clone(self, pop1, pop2):
    return
  
  def _finalize(self):
    pass
  
  def print(self):
    print(self.P)
    print(self.Q)
    print(self.QQ)


class Demography:
  def __init__():
    pass
