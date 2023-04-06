#!/bin/bash
echo "Submitting jobs to generate TennGen data"
export JOB_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm/job-scripts
for i in {0..25} ; do sbatch $JOB_DIR/tenngen-job.sh $1 2 1.1 $i; done