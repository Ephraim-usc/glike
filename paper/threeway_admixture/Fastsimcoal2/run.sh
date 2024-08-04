#!/bin/bash
#SBATCH --ntasks=1
##SBATCH --partition=oneweek
#SBATCH --time=48:00:00
#SBATCH --mem=30000
#SBATCH --array=21-50

name=${SLURM_ARRAY_TASK_ID}

mkdir output_${name}
cd output_${name}


cp ../ThreeWayAdmixture.tpl .
cp ../ThreeWayAdmixture.def .
cp ../ThreeWayAdmixture.est .

if [ ! -f ThreeWayAdmixture_DAFpop0.obs ]; then
echo "simulating data"
fsc28 -t ThreeWayAdmixture.tpl -f ThreeWayAdmixture.def  -n 1 -s0 -d -k 1000000 -q
cp ThreeWayAdmixture/*.obs .
fi


rm ThreeWayAdmixture -r -f
sed -i "$ s/DNA 30000000 1.0e-8 1.0e-8/FREQ 1 0 1.0e-8/g" ThreeWayAdmixture.tpl

for i in `seq 1 20`
do
if [ ! -f ${i}.bestlhoods ]; then
echo "running the ${i}-th estimation"
fsc28 -t ThreeWayAdmixture.tpl -n 100000 -d -e ThreeWayAdmixture.est -M -L 40 -q -c12 -B12
cp ThreeWayAdmixture/ThreeWayAdmixture.bestlhoods ${i}.bestlhoods
rm ThreeWayAdmixture -r -f
fi
done
