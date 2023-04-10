#!/bin/bash
#SBATCH --partition=devel
#SBATCH --job-name=my_conda_job
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=4000


module load miniconda

conda activate qiskitenv3
python 1layer.py