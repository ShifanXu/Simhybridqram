#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --out="/gpfs/gibbs/project/ding/sx86/234010test/slurm_jobs/%a.jobout"
#SBATCH --error="/gpfs/gibbs/project/ding/sx86/234010test/slurm_jobs/%a.joberror"


module load miniconda

conda activate qiskitenv3
python 1layer.py

