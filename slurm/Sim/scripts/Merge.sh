#!/bin/bash

#############################################
# This script is used to generate and submit
# jobs for the ToyModel simulation
#############################################

#############################################
# 1. Set up locations and directories
#############################################
export BASE_PATH=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod # base path
export SLURM_PATH=${BASE_PATH}/slurm/Sim # directory for slurm files
export MACRO_PATH=${BASE_PATH}/macros/sim # directory for the macro
export LOG_PATH=${SLURM_PATH}/logs # directory for log files
export JOB_PATH=${SLURM_PATH}/jobs # directory for job files
export OUTPUT_PATH=${BASE_PATH}/rootfiles # directory for output root files

############ executables
export EXE_NAME=merge_pt_bins # name of the executable


#############################################
# 2. Set up parameters
#############################################

############ arguments
USAGE="Usage: ./Merge.sh <input dir> <name>"
if [ $# -ne 2 ]; then
    echo $USAGE
    exit 1
fi

INPUT_DIR=$1
NAME=$2


# #############################################
# # 4. Clear old logs
# #############################################

# echo "Clearing old logs"
mkdir -p $SLURM_PATH
mkdir -p $LOG_PATH
mkdir -p $JOB_PATH

rm -rf $JOB_PATH/merge_${NAME}*
rm -rf $LOG_PATH/merge_${NAME}*

#############################################
# #############################################
# # 5. Create job files and submit jobs
# #############################################

cd $JOB_PATH
# get nevents per job from the config file

FILENAME=${JOB_PATH}/merge_${NAME}.sh
JOBNAME=merge_${NAME}

echo "#!/bin/bash " >> $FILENAME
echo "#SBATCH -J $JOBNAME			       #The name of the job" >> $FILENAME
echo "#SBATCH -A ACF-UTK0019              # The project account to be charged" >> $FILENAME
echo "#SBATCH --nodes=1                    # Number of nodes" >> $FILENAME
echo "#SBATCH --ntasks-per-node=4         # cpus per node " >> $FILENAME
echo "#SBATCH --partition=condo-cnattras          " >> $FILENAME
# echo "#SBATCH --time=$TIME             # Wall time (days-hh:mm:ss)" >> $FILENAME
echo "#SBATCH --time=00-06:00:00             # Wall time (days-hh:mm:ss)" >> $FILENAME
echo "#SBATCH --error=${LOG_PATH}/${JOBNAME}.e%J	     " >> $FILENAME
echo "#SBATCH --output=${LOG_PATH}/${JOBNAME}.o%J	     " >> $FILENAME
echo "#SBATCH --qos=condo" >> $FILENAME
echo "" >> $FILENAME
echo "" >> $FILENAME
echo "# Load modules and environment" >> $FILENAME
echo "module load anaconda3/2021.05" >> $FILENAME
echo "source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh" >> $FILENAME
echo "conda activate myhepenv" >> $FILENAME
echo "" >> $FILENAME
echo "EXE_DIR=${MACRO_PATH}" >> $FILENAME
echo "EXE=${EXE_NAME}" >> $FILENAME
echo "" >> $FILENAME
echo "INPUT_DIR=${INPUT_DIR}" >> $FILENAME
echo "" >> $FILENAME
echo "echo \"=====================================================\"" >> $FILENAME
echo "echo \"Running \$EXE_DIR/\$EXE with \$INPUT_DIR\"" >> $FILENAME
echo "echo \"=====================================================\"" >> $FILENAME
echo "" >> $FILENAME
echo "cd \$EXE_DIR" >> $FILENAME
echo "./\$EXE \$INPUT_DIR" >> $FILENAME
echo "" >> $FILENAME
echo "echo \"Done\"" >> $FILENAME


# submit job
sbatch_output=$(sbatch $FILENAME)
# get the job id
# Submitted batch job 1335531
job_id=$(echo $sbatch_output | awk '{print $4}')

# echo "Done"