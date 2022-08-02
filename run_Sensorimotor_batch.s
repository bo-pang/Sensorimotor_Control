#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --time=7-00:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=Bo
#SBATCH --mail-type=END
##SBATCH --mail-user=bp1471@nyu.edu
#SBATCH --output=Bo_%A_%a.out
#SBATCH --error=Bo_%A_%a.err

module purge
module load matlab/2017b

cd /scratch/$USER/sensorimotor_rollout
cat HPC_Sensorimotor_batch_$SLURM_ARRAY_TASK_ID.m | matlab -nodisplay