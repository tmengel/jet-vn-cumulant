#!/bin/bash
#SBATCH -J cumulants			       #The name of the job
#SBATCH -A ACF-UTK0019              # The project account to be charged
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=12          # cpus per node 
#SBATCH --partition=condo-cnattras          
#SBATCH --time=0-24:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=cumulants.e%J	     
#SBATCH --output=cumulants.o%J	     
#SBATCH --qos=condo

module load anaconda3/2021.05
source /sw/isaac/applications/anaconda3/2021.05/rhel8_gcc10.2.0/anaconda3-2021.05/etc/profile.d/conda.sh
cd /lustre/isaac/scratch/tmengel/jet-vn-cumulant/analysis
conda activate myhepenv

./calculate_jet_vn /lustre/isaac/scratch/tmengel/jet-vn-cumulant/simulation/root-files/MixedEvents/20-40/200GeV_AlignedMixedEvents_20to40cent_R04.root 0.4 testoutput.root
