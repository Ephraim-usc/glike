import math
import time
import numpy as np
import pandas as pd
import scipy.linalg
from scipy.special import logsumexp
import msprime
np.seterr(divide='ignore')




class Stage: #lineage Markov Process
  def __init__(self, in_pops, out_pops, P, ns):
    self.in_pops = in_pops
    self.out_pops = out_pops
    self.ns = ns
    self.P = P
    self.Q = np.zeros([len(out_pops), len(out_pops)])
    self.QQ = np.zeros([len(out_pops)**2, len(out_pops)**2])
    
  def add_migration(source, dest, rate):
    return
  
  def add_clone(pop1, pop2):
    return
  
  def _finalize():
    


class Demography:
  def __init__():
    pass
