#!/bin/bash
#SBATCH --job-name=progeny
#SBATCH --output=./job_outs/progeny_%j.out
#SBATCH --error=./job_outs/progeny_%j.err
#SBATCH --chdir=/your/parent/dir/
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=1
#SBATCH --time=00-04:00:00
#SBATCH --array=1-4

module load Anaconda3

source /opt/ohpc/pub/libs/easybuild/4.5.0/software/Anaconda3/2022.10/etc/profile.d/conda.sh

conda activate spatial

#Run subset of samples that didn't run the first time due to memory issues
DATA=`sed -n ${SLURM_ARRAY_TASK_ID}p ./Input_Files/samples.txt`
Rscript ./Scripts/04_spatial_progeny_score.R $DATA