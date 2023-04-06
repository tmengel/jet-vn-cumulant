#!/bin/bash
export MACROS_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/Macros
export SLURM_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    echo "Usage: ./merge-all-files.sh <output file> <ptbins> <jet radius>"
    exit 1
fi

echo "output file: $1"
echo "ptbins: $2"
echo "jet radius: $3"

cd $MACROS_DIR
echo "Compiling macros"
make clean
make

echo "Submitting jobs to generate Pythia pp and TennGen data for RHIC"
cd $SLURM_DIR
rm -rf slurm-out/*
rm -rf slurm-err/*


echo "Realigning Pythia pp data for RHIC and merging"
cd $SLURM_DIR/job-scripts 

export MERGED_FILE=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/MixedEvents/20-40/$1 
export INPUT_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/Realigned
sbatch merge-job.sh $INPUT_DIR $2 $3 $MERGED_FILE


echo "Done submitting jobs"

# Path: jet-vn-cumulant/simulation/slurm/generate-pythia-and-tenngen.sh