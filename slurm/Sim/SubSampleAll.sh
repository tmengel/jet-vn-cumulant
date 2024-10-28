#!/bin/bash

input_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/sim
output_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/models
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

sub_dirs=$(ls -d ${input_dir}/*)
rotations=(vn_0 2vn_0 5vn_0 cosine linear logarithmic)
for sub_dir in $sub_dirs ; do
    input_base=${sub_dir}/rotated
    for rotation in ${rotations[@]} ; do
        input=${input_base}/${rotation}
        if [ ! -d ${output_dir}/${rotation} ]; then
            mkdir -p ${output_dir}/${rotation}
        fi
        sub_dir_base=$(basename $sub_dir)
        output=${output_dir}/${rotation}/${sub_dir_base}
        if [ ! -d $output ]; then
            mkdir -p $output
        fi
        cd /lustre/isaac/scratch/tmengel/JetVnCumulantMethod/slurm/Sim/scripts
        ./SubSample.sh $input $output ${sub_dir_base}_${rotation}
    done
done

