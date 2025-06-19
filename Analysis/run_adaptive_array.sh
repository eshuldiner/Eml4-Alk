#!/usr/bin/env bash
#SBATCH --account=mwinslow
#SBATCH -t 48:00:00
#SBATCH --array=1-1
#SBATCH --mem=25g
#SBATCH --requeue
#SBATCH --mail-user=eshuldin
#SBATCH --mail-type=FAIL

#######################
# To run, modify --array flag to match # of samples in project 

# Sample call:
# sbatch run_adaptive_array.sh


#Parameters
parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p $1)

echo $parameters

python3 adaptiveSampling.py $parameters
