#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=30000
#SBATCH --array=1-50

name=${SLURM_ARRAY_TASK_ID}

python3 << END
import os
import msprime
import tsinfer
import tsdate
import glike

demography = glike.threeway_admixture_demography(30, 60, 1e4, 0.4, 0.7, 2000, 20000, 3000, 30000, 10000, 5000)
arg = msprime.sim_ancestry({"O": 1000}, sequence_length = 3e7, recombination_rate = 1e-8, demography = demography, ploidy = 1)
arg = msprime.sim_mutations(arg, rate = 1e-8, discrete_genome = False)
arg.dump("true_${name}.trees")

os.mkdir('tsinfer_${name}')
glike.write_tsinfer_input(arg, 'tsinfer_${name}' + '/tmp')

os.mkdir('relate_${name}')
glike.write_relate_input(arg, 'relate_${name}' + '/tmp')
END


cd tsinfer_${name}
tsinfer infer tmp.samples -p -t 4 --recombination-rate 1e-8
tsdate preprocess tmp.trees preprocessed.trees
tsdate date preprocessed.trees dated.trees 10000 -m 1e-8 --progress
cd ..
cp tsinfer_${name}/dated.trees ./tsdate_${name}.trees
rm -r -f tsinfer_${name}


cd relate_${name}
${PATH_TO_RELATE}/bin/Relate --mode All -m 1e-8 -N 10000 --memory 20 --haps tmp.haps --sample tmp.sample --map tmp.map -o tmp
${PATH_TO_RELATE}/bin/RelateCoalescentRate --mode EstimatePopulationSize -i tmp -o tmp
${PATH_TO_RELATE}/scripts/SampleBranchLengths/SampleBranchLengths.sh -i tmp -o tmp_sampled -m 1e-8 --coal tmp.coal --format a --num_samples 1 --seed 1
${PATH_TO_RELATE}/bin/RelateFileFormats --mode ConvertToTreeSequence -i tmp -o tmp
${PATH_TO_RELATE}/bin/RelateFileFormats --mode ConvertToTreeSequence -i tmp_sampled -o tmp_sampled
cd ..
cp relate_${name}/tmp_sampled.trees ./relate_sampled_${name}.trees
rm -r -f relate_${name}
