#!/bin/bash
base_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/anafiles
sub_dirs=$(ls -d $base_dir/*)
for sub_dir in $sub_dirs; do 
    name=$(basename $sub_dir)
    cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Ana/scripts
    ./Unfold.sh $name ${sub_dir}/MultSubtracted
done