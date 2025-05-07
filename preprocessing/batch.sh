#!/bin/bash

#SBATCH --job-name=preprocessing_ses1
#SBATCH --mail-user=grum90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --time=0-08:00:00
#SBATCH --qos=standard

# add Matlab module
module add MATLAB

# specify task session number (1 or 2)
session=1

# specify preprocessing steps: 
# 4=realignment; 5=coregistration; 1=segmentation; 6=normalization; 8=CSFnoise correction; 9=smoothing
analysis_steps="[4, 5, 1, 8]"
# initial prefix for the data (empty if raw)
prefix=""

# batch input as subject number for parallelization
subject=${SLURM_ARRAY_TASK_ID}

# run function
matlab -r "preprocessing_hpc(${subject}, ${session}, ${analysis_steps}, '${prefix}')" $MCR_HOME 

