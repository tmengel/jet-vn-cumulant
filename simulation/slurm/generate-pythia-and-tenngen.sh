#!/bin/bash
export MACROS_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/Macros
export SLURM_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: ./generate-pythia-and-tenngen.sh <energy> <number of events>"
    exit 1
fi

echo "Input energy (GeV): $1"
echo "Number of events: $2"

cd $MACROS_DIR
echo "Compiling macros"
make clean
make

echo "Submitting jobs to generate Pythia pp and TennGen data for RHIC"
cd $SLURM_DIR
rm -rf slurm-out/*
rm -rf slurm-err/*


echo "Submitting jobs to generate Pythia pp data for RHIC"
cd shells 
./submit-pythia-slurm-jobs.sh $1 $2
echo "Submitting jobs to generate TennGen data for RHIC"
./submit-tenngen-slurm-jobs.sh $2 

echo "Done submitting jobs"

# Path: jet-vn-cumulant/simulation/slurm/generate-pythia-and-tenngen.sh