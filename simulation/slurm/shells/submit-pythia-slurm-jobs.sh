#!/bin/bash
echo "Submitting jobs to generate Pythia pp data"
export JOB_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm/job-scripts
for i in {0..25}; do sbatch $JOB_DIR/pythia-job.sh $1 $2 $i; done
