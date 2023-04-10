#!/bin/bash
#SBATCH --job-name=my_job1
#SBATCH --partition=day
#SBATCH --time=4:00:00
#SBATCH --mem=160G
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1S
#SBATCH --cpus-per-task=48
#SBATCH --out="/gpfs/gibbs/project/ding/sx86/234010test/a.jobout"
#SBATCH --error="/gpfs/gibbs/project/ding/sx86/234010test/a.joberror"


module load miniconda

conda activate qiskitenv3
python reduced_2layer.py