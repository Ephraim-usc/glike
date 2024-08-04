#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=oneweek
#SBATCH --time=96:00:00
#SBATCH --mem=30000
#SBATCH --array=1-100

files=(*.trees)
file=${files[$((SLURM_ARRAY_TASK_ID-1))]}
name=${file%%.*}


python3 << END
from glike import *

print("input file: ${name}.trees", flush = True)
arg = tskit.load("${name}.trees")
trees = [arg.at(pos).copy() for pos in range(300000, 30000000, 300000)]

def fun(t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix):
  tmp = ["admix"] * 1000
  samples = {i:pop for i, pop in enumerate(tmp)}
  demo = glike.threeway_admixture_demo(t1, t2, t3, t4, r1, r2, N_afr, N_eur, N_asia, N_admix, N_ooa, N_anc, gr_eur, gr_asia, gr_admix)
  return glike.glike_trees(trees, demo, samples = samples, prune = 0.5)

x0 = {"t1":5, "t2":200, "t3":500, 
      "t4":2000, "r1":0.33, "r2":0.33, 
      "N_afr":10000, "N_eur":10000, "N_asia":10000, "N_admix":10000, "N_ooa":10000, "N_anc":10000,
      "gr_eur":0.01, "gr_asia":0.01, "gr_admix":0.01}
bounds = [(1,"t2"),("t1","t3"),("t2","t4"),("t3",100000),(0.0, "1-r2"),(0.0, "1-r1"),(100,1e6),(100,1e6),(100,1e6),(100,1e6),(100,1e6),(100,1e6),(0.0001,1),(0.0001,1),(0.0001,1)]
x, logp = glike.maximize(fun, x0, bounds = bounds)

with open("${name}.txt", "at") as out_file:
  out_file.write("${name}" + "\t" + "american_admixture_demo" + "\t" + str(x) + "\t" + str(logp) + "\n")
END
