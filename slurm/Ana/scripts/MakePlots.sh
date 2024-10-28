#!/bin/bash
plot_exe=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/plots/PlotFlow.py
output_path=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/plots/result_May_16_2024_5to30toInfGeV_50MillJetEvents
# Rotations=(ConstantA ConstantB ConstantC Linear Cosine Logarithmic)
Rotations=(ConstantBLarge)
file_dir=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod/rootfiles/ConstantBLarge
files=()
for i in "${Rotations[@]}"
do
    files+=(${file_dir}/${i}/flow_results.root)
done


$plot_exe -o $output_path -f ${files[@]} -r ${Rotations[@]}

