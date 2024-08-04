#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=oneweek
#SBATCH --time=168:00:00
#SBATCH --mem=20000
#SBATCH --array=21-50

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import numpy as np
import msprime
import tskit
import glike
import math
import pickle
from time import time as tt

model_demography = glike.ancient_europe_demography
model_demo = glike.ancient_europe_demo
truth = model_demography.__defaults__
print(f"truth = {truth}\n", flush = True)

samples = [msprime.SampleSet(20, population = "ana", time = 260),
           msprime.SampleSet(20, population = "neo", time = 180),
           msprime.SampleSet(20, population = "whg", time = 250),
           msprime.SampleSet(100, population = "bronze", time = 0), 
           msprime.SampleSet(20, population = "yam", time = 160),
           msprime.SampleSet(20, population = "ehg", time = 250),
           msprime.SampleSet(20, population = "chg", time = 300)]

demography = model_demography(*truth)
arg = msprime.sim_ancestry(samples, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
trees = [arg.at(pos).copy() for pos in range(300000, 30000000, 300000)]


def fun(t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr):
  tmp = ["ana"] * 20 + ["neo"] * 20 + ["whg"] * 20 + ["bronze"] * 100 + ["yam"] * 20 + ["ehg"] * 20 + ["chg"] * 20
  samples = {i:pop for i, pop in enumerate(tmp)}
  demo = glike.ancient_europe_demo(t1, t2, t3, t4, t5, t6, r1, r2, r3, N_ana, N_neo, N_whg, N_bronze, N_yam, N_ehg, N_chg, N_ne, N_wa, N_ooa, gr)
  return glike.glike_trees(trees, demo, samples = samples, prune = 0.5)

x0 = {"t1":120, "t2":170, "t3":190, "t4":500, "t5":600, "t6":1300, 
      "r1":0.5, "r2":0.5, "r3":0.5, 
      "N_ana":10000, "N_neo":10000, "N_whg":10000, "N_bronze":10000, "N_yam":10000, "N_ehg":10000, "N_chg":10000, "N_ne":10000, "N_wa":10000, "N_ooa":10000, 
      "gr":0.01}
bounds = [(1,"t2"),("t1","t3"),("t2","t4"),("t3", "t5"),("t4","t6"),("t5",10000),
          (0.01, 0.99),(0.01,0.99),(0.01,0.99),
          (100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),
          (0.0001,1.0)]
x, logp = glike.maximize(fun, x0, bounds = bounds)

with open("${name}.txt", "at") as out_file:
  out_file.write("${name}" + "\t" + "ancient_europe_demo" + "\t" + str(x) + "\t" + str(logp) + "\n")
END
