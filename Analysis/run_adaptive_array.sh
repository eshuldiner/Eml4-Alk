#!/usr/bin/env bash
#SBATCH --array=1-2

#######################
# To run, modify --array flag to match # of samples in project 

# Sample call:
# sbatch run_adaptive_array.sh

#Parameters
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p $1)

echo $parameters

python3 adaptiveSampling.py $parameters
