#!/bin/bash
#SBATCH --job-name=argo_test
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/argo_test.out
#SBATCH --error=logs/argo_test.err

module load jaspy

echo "Starting Argo test job"
python scripts/argo_mean_temp.py
echo "Job finished"