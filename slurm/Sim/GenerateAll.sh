#!/bin/bash
iptBin=$1
ptbins=(5to10 10to15 15to20 20to25 25to30 30to35 35to40 40to50 50to60 60to80)
nJobs=100
# done : 0 1 2 3 4 5 6 7
# to do : 8 9
config_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/configs
config_file=$(ls ${config_dir}/*${ptbins[iptBin]}*.txt)
echo "Generating jobs for $config_file"
cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
./Generate.sh $config_file $nJobs


