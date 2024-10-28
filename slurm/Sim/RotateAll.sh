#!/bin/bash

config_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/configs/rotation
for fileidx in {0..100}; do 
    # echo "Rotation jobs for $fileidx"
    cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
    ./Rotate.sh $config_dir $fileidx
done


