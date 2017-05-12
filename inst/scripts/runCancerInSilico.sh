#!/bin/bash -l

#SBATCH
#SBATCH --job-name=CancerInSilico
#SBATCH --time=5:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=rcao5@jhu.edu

if [ "$1" != "" ]; then
    jobNum=$1
else
    echo missing argument
    exit 1
fi

if [ "$SLURM_ARRAY_TASK_ID" != "" ]; then
    arrayNum=$SLURM_ARRAY_TASK_ID
else
    arrayNum=1
fi

time Rscript ~/CancerInSilico/inst/scripts/runCancerInSilico.R $arrayNum $jobNum
