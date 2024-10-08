#!/bin/bash
#SBATCH --ntasks=1
##SBATCH --partition=oneweek
#SBATCH --time=48:00:00
#SBATCH --mem=30000
#SBATCH --array=1-50

name=${SLURM_ARRAY_TASK_ID}

mkdir output_${name}
cd output_${name}


cp ../AncientEurope.tpl .
cp ../AncientEurope.def .
cp ../AncientEurope.est .

if [ ! -f AncientEurope_DAFpop0.obs ]; then
echo "simulating data"
fsc28 -t AncientEurope.tpl -f AncientEurope.def  -n 1 -s0 -d -k 10000000 -q --multiSFS
cp AncientEurope/*.obs .
fi


rm AncientEurope -r -f
sed -i "$ s/DNA 30000000 1.0e-8 1.0e-8/FREQ 1 0 1.0e-8/g" AncientEurope.tpl

for i in `seq 1 20`
do
if [ ! -f ${i}.bestlhoods ]; then
echo "running the ${i}-th estimation"
fsc28 -t AncientEurope.tpl -n 100000 -d -e AncientEurope.est -M -L 40 -q --multiSFS -c12 -B12
cp AncientEurope/AncientEurope.bestlhoods ${i}.bestlhoods
rm AncientEurope -r -f
fi
done
