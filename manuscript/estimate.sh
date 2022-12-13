#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=20000
#SBATCH --array=1-60

files=(*.trees)
file=${files[$SLURM_ARRAY_TASK_ID]}
name=${file%%.*}


python3 << END
from glike import *
N = 1000
l = 30000000

arg = tskit.load("${name}.trees")
trees = [arg.at(pos).copy() for pos in range(300000, l, 300000)]

names = ["t1", "t2", "r", "N", "N_a", "N_b", "N_c", "m_ab"]
values = [20, 50000.0, 0.5, 5000, 20000, 20000, 5000, 0, 0]
limits = [(1,100),(1e4, 2e5),(0.001,0.999),(100,100000),(100,100000),(100,100000),(100,100000),(0,0)]
fixed = ["m_ab"]
search = Search(names, values, limits, fixed)
x, logp = estimate(trees, twoway_admixture_demo, search)
with open("${name}.txt", "at") as out_file:
  out_file.write("${name}" + "\t" + "twoway_admixture_demo" + "\t" + str(x) + "\t" + str(logp) + "\n")

names = ["t1", "t2", "t3", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e", "m_ab", "m_cd"]
values = [20, 80, 50000.0, 0.5, 0.5, 5000, 20000, 5000, 20000, 20000, 5000, 0, 0]
limits = [(1,"t2"),("t1",100),(1e4, 2e5),(0.001,0.999),(0.001,0.999),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(0,0),(0,0)]
fixed = ["m_ab", "m_cd"]
search = Search(names, values, limits, fixed)
x, logp = estimate(trees, threeway_admixture_demo, search)
with open("${name}.txt", "at") as out_file:
  out_file.write("${name}" + "\t" + "threeway_admixture_demo" + "\t" + str(x) + "\t" + str(logp) + "\n")

END
