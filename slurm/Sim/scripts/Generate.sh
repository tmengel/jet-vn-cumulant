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
export QA_PATH=${SLURM_PATH}/QA # directory for QA files

############ executables
export EXE_NAME=generate # name of the executable


#############################################
# 2. Set up parameters
#############################################

############ arguments
USAGE="Usage: ./Generate.sh <config file> <nJobs>"
if [ $# -ne 2 ]; then
    echo $USAGE
    exit 1
fi

CONFIG_FILE=$1
N_JOBS=$2

# make sure the config file exists
if [ ! -f $CONFIG_FILE ]; then
    echo "Config file $CONFIG_FILE does not exist in $CONFIG_PATH"
    exit 1
fi


CONFIG_BASENAME=$(basename $CONFIG_FILE)
# remove the extension
CONFIG_DIR=$(dirname $CONFIG_FILE)
cd $CONFIG_DIR
CONFIG_ABS_PATH=$(pwd)/$CONFIG_BASENAME
CONFIG_BASENAME=${CONFIG_BASENAME%.*}

# echo "Config file: $CONFIG_FILE"
echo "Config path: $CONFIG_ABS_PATH"
echo "Config basename: $CONFIG_BASENAME"
NEVENTS_PER_JOB=$(grep "Main::NumEvents" $CONFIG_ABS_PATH | awk -F "==" '{print $2}')
NEVENTS_PER_JOB=$(echo $NEVENTS_PER_JOB | xargs)
echo "Number of events per job: $NEVENTS_PER_JOB"

############################################
# 3. Compile the executable
############################################

# echo "Compiling executable"
# cd $EXE_PATH
# make clean 
# make

# #############################################
# # 4. Clear old logs
# #############################################

# echo "Clearing old logs"
mkdir -p $SLURM_PATH
mkdir -p $LOG_PATH
mkdir -p $JOB_PATH

rm -rf $JOB_PATH/${CONFIG_BASENAME}*.*
rm -rf $LOG_PATH/${CONFIG_BASENAME}*.*

#############################################
# 5. Clear old root files
#############################################

export OUTPUT_PATH=${OUTPUT_PATH}/sim
mkdir -p $OUTPUT_PATH
OUTPUT_PATH=${OUTPUT_PATH}/${CONFIG_BASENAME}

# remove old files if exists
if [ -d $OUTPUT_PATH ]; then
    rm -rf $OUTPUT_PATH
fi
mkdir -p $OUTPUT_PATH

BKGD_OUTPUT_PATH=${OUTPUT_PATH}/bkgd
mkdir -p $BKGD_OUTPUT_PATH
SIG_OUTPUT_PATH=${OUTPUT_PATH}/sig
mkdir -p $SIG_OUTPUT_PATH
#############################################
# 6. Setup QA output
#############################################
# make sure the QA directory exists
if [ ! -d $QA_PATH ]; then
    mkdir -p $QA_PATH
fi

# delete old QA files 
rm -rf $QA_PATH/${CONFIG_BASENAME}*.txt

# create a QA file 
QA_FILE=$QA_PATH/${CONFIG_BASENAME}_QA.txt

echo "QA file: $QA_FILE"
echo "Config file: $CONFIG_ABS_PATH" > $QA_FILE
echo "Number of jobs: $N_JOBS" >> $QA_FILE
echo "Number of events per job: $NEVENTS_PER_JOB" >> $QA_FILE
echo "Output path: $OUTPUT_PATH" >> $QA_FILE
echo "" >> $QA_FILE

#############################################
# #############################################
# # 6. Create job files and submit jobs
# #############################################

cd $JOB_PATH
BKGD_OUTFILE_BASE=${BKGD_OUTPUT_PATH}/${CONFIG_BASENAME}_job_
SIG_OUTFILE_BASE=${SIG_OUTPUT_PATH}/${CONFIG_BASENAME}_job_
# get nevents per job from the config file
for i in $(seq 0 $(($N_JOBS-1))); do 
    
    # write job file
    FILENAME=${JOB_PATH}/${CONFIG_BASENAME}_job_${i}.sh
    JOBNAME=${CONFIG_BASENAME}_${i}
    BKGD_OUTFILE=${BKGD_OUTFILE_BASE}${i}.root
    SIG_OUTFILE=${SIG_OUTFILE_BASE}${i}.root

    echo "#!/bin/bash " >> $FILENAME
    echo "#SBATCH -J $JOBNAME			       #The name of the job" >> $FILENAME
    echo "#SBATCH -A ACF-UTK0019              # The project account to be charged" >> $FILENAME
    echo "#SBATCH --nodes=1                    # Number of nodes" >> $FILENAME
    echo "#SBATCH --ntasks-per-node=1         # cpus per node " >> $FILENAME
    echo "#SBATCH --partition=condo-cnattras          " >> $FILENAME
    # echo "#SBATCH --time=$TIME             # Wall time (days-hh:mm:ss)" >> $FILENAME
    echo "#SBATCH --time=00-04:00:00             # Wall time (days-hh:mm:ss)" >> $FILENAME
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
    echo "CONFIG_FILE=$CONFIG_ABS_PATH" >> $FILENAME
    echo "OUTPUT_BKG_FILE=$BKGD_OUTFILE" >> $FILENAME
    echo "OUTPUT_SIG_FILE=$SIG_OUTFILE" >> $FILENAME
    echo "" >> $FILENAME
    echo "echo \"=====================================================\"" >> $FILENAME
    echo "echo \"Running \$EXE_DIR/\$EXE with \$CONFIG_FILE\"" >> $FILENAME
    echo "echo \"Background output: \$OUTPUT_BKG_FILE\"" >> $FILENAME
    echo "echo \"Signal output: \$OUTPUT_SIG_FILE\"" >> $FILENAME
    echo "echo \"=====================================================\"" >> $FILENAME
    echo "# wait for $i seconds" >> $FILENAME
    echo "sleep $i" >> $FILENAME
    echo "" >> $FILENAME
    echo "cd \$EXE_DIR" >> $FILENAME
    echo "./\$EXE \$CONFIG_FILE \$OUTPUT_SIG_FILE \$OUTPUT_BKG_FILE" >> $FILENAME
    echo "" >> $FILENAME
    echo "echo \"Done\"" >> $FILENAME


    # submit job
    sbatch_output=$(sbatch $FILENAME)
    # get the job id
    # Submitted batch job 1335531
    job_id=$(echo $sbatch_output | awk '{print $4}')

    # write to QA file
    echo "Job$i: $JOBNAME, ${FILENAME}, ${BKGD_OUTFILE}, ${SIG_OUTFILE}, ${LOG_PATH}/${JOBNAME}.e${job_id}, ${LOG_PATH}/${JOBNAME}.e${job_id}, $job_id" >> $QA_FILE

done


echo "Done"