#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=30000
#SBATCH --array=1-20

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import os
from glike import *
import tsinfer
import tsdate

N = 1000
l = 30000000

demography = threeway_admixture_demography(30, 60, 1e5, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000, 3e-4, 1e-4)
#demography = twoway_admixture_demography(40, 1e5, 0.3, 3000, 20000, 10000, 8000, 0)
arg = msprime.sim_ancestry({"O": N}, sequence_length = l, recombination_rate = 1e-10, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("true_${name}.trees")

os.mkdir('tsinfer_${name}')
write_tsinfer_input(arg, 'tsinfer_${name}' + '/rec')

os.mkdir('relate_${name}')
write_relate_input(arg, 'relate_${name}' + '/rec')
END


cd tsinfer_${name}
~/bin/tsinfer infer rec.samples -p -t 4 --recombination-rate 1e-10 
~/bin/tsdate preprocess rec.trees preprocessed.trees
~/bin/tsdate date preprocessed.trees dated.trees 14834 -m 1e-8 --progress
cd ..
cp tsinfer_${name}/dated.trees ./tsdate_${name}.trees
rm -r -f tsinfer_${name}


cd relate_${name}
~/bin/Relate --mode All -m 1e-8 -N 34878 --memory 20 --haps rec.haps --sample rec.sample --map rec.map -o rec
~/bin/RelateFileFormats --mode ConvertToTreeSequence -i rec -o rec
cd ..
cp relate_${name}/rec.trees ./relate_${name}.trees
rm -r -f relate_${name}
