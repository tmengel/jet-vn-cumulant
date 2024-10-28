#!/bin/bash
ptBins=(5to10 10to15 15to20 20to25 25to30 30to35 35to40 40to50 50to60 60to80)
nJobs=100
cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
./GenerateAll.sh 0
./GenerateAll.sh 1

root_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/sim
for ptbins in 2 3 4 5 6 7 8 9 ; do 
    previous_ptbins=$((ptbins-2))

    dir=${root_dir}/toy_model_config_${ptBins[previous_ptbins]}_GeV/bkgd
    while [ ! -d $dir ]; do
        echo "Waiting for the previous job to finish"
        sleep 300
    done

    nfiles=$(ls -l ${dir}/*.root | wc -l)
    echo "nfiles = $nfiles"
    while [ $nfiles -lt $nJobs ]; do
        echo "Waiting for the previous job to finish"
        sleep 300
        nfiles=$(ls -l ${dir}/*.root | wc -l)
    done
    echo "Generating jobs for $config_file"
    cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
    ./GenerateAll.sh $config_file $nJobs
done


