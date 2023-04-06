#!/bin/bash
echo "Submitting jobs to align"
export JOB_DIR=/lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/slurm/job-scripts

export PYTHIA_FILE_BASE=$1/200GeV_PP_ptbin
export TENNGEN_FILE_BASE=$2/200GeV_AuAu_ptbin

for i in {0..25} ; do
    export PYTHIA_FILE=$PYTHIA_FILE_BASE$i.root
    export TENNGEN_FILE=$TENNGEN_FILE_BASE$i.root
    # echo $PYTHIA_FILE
    # echo $TENNGEN_FILE
    # echo $3
    # echo $4
    sbatch $JOB_DIR/realign-job.sh $PYTHIA_FILE $TENNGEN_FILE $3 $4;
done