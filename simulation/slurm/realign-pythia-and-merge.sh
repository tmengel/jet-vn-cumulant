#!/bin/bash
export MACROS_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/Macros
export SLURM_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm

# if [ $# -eq 0 ]
#   then
#     echo "No arguments supplied"
#     echo "Usage: ./realign-pythia-and-merge.sh <energy> <number of events>"
#     exit 1
# fi

# echo "Input energy (GeV): $1"
# echo "Number of events: $2"

cd $MACROS_DIR
echo "Compiling macros"
make clean
make

echo "Submitting jobs to generate Pythia pp and TennGen data for RHIC"
cd $SLURM_DIR
rm -rf slurm-out/*
rm -rf slurm-err/*


echo "Realigning Pythia pp data for RHIC and merging"
cd shells 

export PYTHIA_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/PP
export TENNGEN_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/AuAu/20-40
export MERGED_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/Realigned
./submit-realign-slurm-jobs.sh $PYTHIA_DIR $TENNGEN_DIR $MERGED_DIR 0.4 100

echo "Done submitting jobs"

# Path: jet-vn-cumulant/simulation/slurm/generate-pythia-and-tenngen.sh