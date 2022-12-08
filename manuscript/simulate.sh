#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=30000
#SBATCH --array=1-10

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import os
from glike import *
import tsinfer
import tsdate

N = 1000
l = 30000000

demography = threeway_admixture_demography(30, 60, 1e5, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000, 0, 0)
arg = msprime.sim_ancestry({"O": N}, sequence_length = l, recombination_rate = 1e-10, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("true_${name}.trees")

trees = [arg.at(pos).copy() for pos in range(0, l, 1000000)]
names = ["t1", "t2", "t3", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e", "m_ab", "m_cd"]
values = [10, 30, 5e4, 0.5, 0.5, 10000, 10000, 10000, 10000, 10000, 10000, 0, 0]
limits = [(0,"t2"),("t1",100),(1e4, 2e5),(0,1),(0,1),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(0,0),(0,0)]
fixed = ["m_ab", "m_cd"]

searchspace = Searchspace(names, values, limits, fixed)
estimate(trees, threeway_admixture_demo, searchspace)

os.mkdir('tsinfer_${name}')
write_tsinfer_input(arg, 'tsinfer_${name}' + '/rec')

os.mkdir('relate_${name}')
write_relate_input(arg, 'relate_${name}' + '/rec')
END


cd tsinfer_${name}
~/bin/tsinfer infer rec.samples -p -t 4 --recombination-rate 1e-10
~/bin/tsdate preprocess rec.trees preprocessed.trees
~/bin/tsdate date preprocessed.trees dated.trees 14834 -m 1e-8 --progress


python3 << END
from glike import *
N = 1000
l = 30000000

arg_tsdate = tskit.load("dated.trees")
trees_tsdate = [arg_tsdate.at(pos).copy() for pos in range(1000000, l, 1000000)]

names = ["t1", "t2", "t3", "r1", "r2", "N", "N_a", "N_b", "N_c", "N_d", "N_e", "m_ab", "m_cd"]
values = [10, 30, 5e4, 0.5, 0.5, 10000, 10000, 10000, 10000, 10000, 10000, 0, 0]
limits = [(0,"t2"),("t1",100),(1e4, 2e5),(0,1),(0,1),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(100,100000),(0,0),(0,0)]
fixed = ["m_ab", "m_cd"]

searchspace = Searchspace(names, values, limits, fixed)
estimate(trees_tsdate, threeway_admixture_demo, searchspace)
END
cd ..

cd relate_${name}
~/bin/Relate --mode All -m 1e-8 -N 34878 --memory 20 --haps rec.haps --sample rec.sample --map rec.map -o rec
~/bin/RelateFileFormats --mode ConvertToTreeSequence -i rec -o rec
cd ..











