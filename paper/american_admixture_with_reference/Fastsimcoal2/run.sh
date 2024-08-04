#!/bin/bash
#SBATCH --ntasks=1
##SBATCH --partition=oneweek
#SBATCH --time=48:00:00
#SBATCH --mem=30000
#SBATCH --array=21-50

name=${SLURM_ARRAY_TASK_ID}

mkdir output_${name}
cd output_${name}


cp ../AmericanAdmixture.tpl .
cp ../AmericanAdmixture.def .
cp ../AmericanAdmixture.est .

if [ ! -f AmericanAdmixture_DAFpop0.obs ]; then
echo "simulating data"
fsc28 -t AmericanAdmixture.tpl -f AmericanAdmixture.def  -n 1 -s0 -d -k 1000000 -q --multiSFS
cp AmericanAdmixture/*.obs .
fi

rm AmericanAdmixture -r -f
sed -i "$ s/DNA 30000000 1.0e-8 1.0e-8/FREQ 1 0 1.0e-8/g" AmericanAdmixture.tpl

for i in `seq 1 20`
do
if [ ! -f ${i}.bestlhoods ]; then
echo "running the ${i}-th estimation"
fsc28 -t AmericanAdmixture.tpl -n 100000 -d -e AmericanAdmixture.est -M -L 40 -q --multiSFS -c12 -B12
cp AmericanAdmixture/AmericanAdmixture.bestlhoods ${i}.bestlhoods
rm AmericanAdmixture -r -f
fi
done
