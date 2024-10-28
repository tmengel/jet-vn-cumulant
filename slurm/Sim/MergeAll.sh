#!/bin/bash

input_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/models
sub_dirs=$(ls -d ${input_dir}/*)
for sub_dir in $sub_dirs ; do
    name=$(basename $sub_dir)
    cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
    ./Merge.sh $sub_dir $name
done

