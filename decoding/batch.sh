#!/bin/bash

#SBATCH --job-name=decoding_ses1
#SBATCH --mail-user=grum90@zedat.fu-berlin.de
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --time=1-12:00:00 
#SBATCH --qos=standard

# add Matlab module
module add MATLAB

# specify task session number (1 or 2)
session=1

# specify preprocessing steps: 
# 1=spec conditions; 2=FIR; 3=GLM; 4=decoding; 5=average; 6=normalize; 7=smooth
analysis_steps="[1 2 4 5 6 7]"
# initial prefix for the data 
prefix="r"

# batch input as subject number for parallelization
subject=${SLURM_ARRAY_TASK_ID}

matlab -r "decoding_hpc(${subject}, ${session}, ${analysis_steps}, '${prefix}')" $MCR_HOME 



