#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=20000
#SBATCH --array=1-50

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import os
import math
import tsinfer
import tsdate
import glike

N_eur = 1032 * math.exp(0.0038 * 920)
N_asia = 554 * math.exp(0.0048 * 920)
N_admix = 30000 * math.exp(0.05 * 12)

demography = glike.american_admixture_demography(12, 920, 2040, 5920, 0.17, 0.33, 14474, N_eur, N_asia, N_admix, 1861, 7310, 0.0038, 0.0048, 0.05)
arg = msprime.sim_ancestry({"admix":1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("true_${name}.trees")

os.mkdir('tsinfer_${name}')
glike.write_tsinfer_input(arg, 'tsinfer_${name}' + '/rec')
END


cd tsinfer_${name}
~/bin/tsinfer infer rec.samples -p -t 4 --recombination-rate 1e-10 
~/bin/tsdate preprocess rec.trees preprocessed.trees
~/bin/tsdate date preprocessed.trees dated.trees 10000 -m 1e-8 --progress
cd ..
cp tsinfer_${name}/dated.trees ./tsdate_${name}.trees
rm -r -f tsinfer_${name}
