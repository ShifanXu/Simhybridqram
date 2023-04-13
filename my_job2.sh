#!/bin/bash
#SBATCH --job-name=my_job2
#SBATCH --partition=day
#SBATCH --time=04:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --out="/gpfs/gibbs/project/ding/sx86/234010test/a.jobout"
#SBATCH --error="/gpfs/gibbs/project/ding/sx86/234010test/a.joberror"


module load miniconda

conda activate qiskitenv3
python reduced_12layer.py