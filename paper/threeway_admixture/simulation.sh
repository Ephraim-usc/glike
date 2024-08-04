#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=30000
#SBATCH --array=1-50

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import os
import glike
import tsinfer
import tsdate

demography = glike.threeway_admixture_demography(30, 60, 1e5, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000)
arg = msprime.sim_ancestry({"O": 1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("true_${name}.trees")

os.mkdir('tsinfer_${name}')
glike.write_tsinfer_input(arg, 'tsinfer_${name}' + '/rec')

os.mkdir('relate_${name}')
glike.write_relate_input(arg, 'relate_${name}' + '/rec')
END


cd tsinfer_${name}
~/bin/tsinfer infer rec.samples -p -t 4 --recombination-rate 1e-8
~/bin/tsdate preprocess rec.trees preprocessed.trees
~/bin/tsdate date preprocessed.trees dated.trees 10000 -m 1e-8 --progress
cd ..
cp tsinfer_${name}/dated.trees ./tsdate_${name}.trees
rm -r -f tsinfer_${name}


cd relate_${name}
~/bin/Relate --mode All -m 1e-8 -N 40000 --memory 20 --haps rec.haps --sample rec.sample --map rec.map -o rec
~/bin/RelateFileFormats --mode ConvertToTreeSequence -i rec -o rec
cd ..
cp relate_${name}/rec.trees ./relate_${name}.trees
rm -r -f relate_${name}
