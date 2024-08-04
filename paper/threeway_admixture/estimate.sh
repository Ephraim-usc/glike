#!/bin/bash
#SBATCH --begin=now+0hour
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=30000
#SBATCH --array=1-150

files=(*.trees)
file=${files[$(expr $SLURM_ARRAY_TASK_ID - 1)]}
name=${file%%.*}


python3 << END
import tskit
import glike

arg = tskit.load("${name}.trees")
trees = [arg.at(pos).copy() for pos in range(300000, 30000000, 300000)]

def fun(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e):
  demo = glike.threeway_admixture_demo(t1, t2, t3, r1, r2, N, N_a, N_b, N_c, N_d, N_e)
  return glike.glike_trees(trees, demo, prune = 0.5)

x0 = {"t1":10, "t2":20, "t3":5000.0, "r1":0.5, "r2":0.5, "N":10000, "N_a":10000, "N_b":10000, "N_c":10000, "N_d":10000, "N_e":10000}
bounds = [(1, "t2"), ("t1", "t3"), ("t2", 1e5), (0.0,1.0), (0.0,1.0), (100,1000000), (100,1000000), (100,1000000), (100,1000000), (100,1000000), (100,1000000)]
x, logp = glike.maximize(fun, x0, bounds = bounds)

with open("${name}.txt", "at") as out_file:
  out_file.write("${name}" + "\t" + "threeway_admixture_demo" + "\t" + str(x) + "\t" + str(logp) + "\n")
END
