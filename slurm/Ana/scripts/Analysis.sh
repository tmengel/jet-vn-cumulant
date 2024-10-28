#!/bin/bash

#############################################
# 1. Set up locations and directories
#############################################
export BASE_PATH=/lustre/isaac/scratch/tmengel/JetVnCumulantMethod # base path
export SLURM_PATH=${BASE_PATH}/slurm/Ana # directory for slurm files
export EXE_PATH=${BASE_PATH}/macros/ana # directory for the macro
export LOG_PATH=${SLURM_PATH}/logs # directory for log files
export JOB_PATH=${SLURM_PATH}/jobs # directory for job files
export OUTPUT_PATH=${BASE_PATH}/rootfiles/anafiles # directory for output root files
export SRC_PATH=${BASE_PATH}/jet-vn-cumulant/analysis/src # directory for the source code

############ executables
export EXE_NAME=analysis # name of the executable


ROTATION=$1
INPUT_DIR=$2


if [ -d $OUTPUT_PATH/$ROTATION ]; then
    rm -rf $OUTPUT_PATH/$ROTATION
fi
mkdir -p $OUTPUT_PATH/$ROTATION
export OUTPUT_DIR=$OUTPUT_PATH/$ROTATION
# make library, exit if error
# pipe output to /dev/null to avoid printing to screen
cd $SRC_PATH
make > /dev/null
# exit if error
if [ $? -ne 0 ]; then
    echo "Error: make failed"
    exit 1
fi

cd $EXE_PATH
make > /dev/null

# echo "Clearing old logs"
cd $JOB_PATH
mkdir -p $SLURM_PATH
mkdir -p $LOG_PATH
mkdir -p $JOB_PATH

rm -rf ${JOB_PATH}/${ROTATION}*.*
rm -rf ${LOG_PATH}/${ROTATION}_ana*.*

FILENAME=${JOB_PATH}/${ROTATION}_ana.sh
JOBNAME=${ROTATION}_analysis

echo "#!/bin/bash " >> $FILENAME
echo "#SBATCH -J $JOBNAME			       #The name of the job" >> $FILENAME
echo "#SBATCH -A ACF-UTK0019              # The project account to be charged" >> $FILENAME
echo "#SBATCH --nodes=1                    # Number of nodes" >> $FILENAME
echo "#SBATCH --ntasks-per-node=10         # cpus per node " >> $FILENAME
echo "#SBATCH --partition=condo-cnattras          " >> $FILENAME
echo "#SBATCH --time=00-05:00:00             # Wall time (days-hh:mm:ss)" >> $FILENAME
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
echo "EXE_DIR=$EXE_PATH" >> $FILENAME
echo "EXE=$EXE_NAME" >> $FILENAME
echo "" >> $FILENAME
echo "" >> $FILENAME
echo "OUTPUT_DIR=$OUTPUT_DIR" >> $FILENAME
echo "INPUT_DIR=$INPUT_DIR" >> $FILENAME
echo "echo \"=====================================================\"" >> $FILENAME
echo "echo \"Running \$EXE_DIR/\$EXE \"" >> $FILENAME
echo "echo \"Output directory: \$OUTPUT_DIR\"" >> $FILENAME
echo "echo \"Input dir: \$INPUT_DIR\"" >> $FILENAME 
echo "echo \"=====================================================\"" >> $FILENAME
echo "" >> $FILENAME
echo "cd \$EXE_DIR" >> $FILENAME
echo "./\$EXE \$INPUT_DIR \$OUTPUT_DIR" >> $FILENAME
echo "" >> $FILENAME
echo "" >> $FILENAME
echo "echo \"=====================================================\"" >> $FILENAME
echo "echo \"=====================================================\"" >> $FILENAME
echo "echo \"Done\"" >> $FILENAME

sbatch_output=$(sbatch $FILENAME)

# echo $sbatch_output
echo "Job submitted: $ROTATION"
