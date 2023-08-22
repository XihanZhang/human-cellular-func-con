#!/bin/bash
#SBATCH --job-name=SpinVIS
#SBATCH -p psych_week
#SBATCH --out="slurm-%j.out"
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2
#SBATCH --mem-per-cpu=50G
#SBATCH --mail-type=ALL

# Set up the environment
module load miniconda
conda activate spintest-environment

# Run the python script
python /gpfs/milgram/project/holmes/xz555/gradient_shift/scripts/08a_SpintestCellTypeNull_VIS.py