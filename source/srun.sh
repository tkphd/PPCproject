#!/bin/bash
#
# srun command for CSCI-6360 group project.
# USAGE: /full/path/to/./q_GG.out [--help] [--init dimension [outfile]] [infile [outfile] steps [increment]]
#
#SBATCH --mail-type=END
#SBATCH --mail-user=kellet@rpi.edu
#SBATCH -D /gpfs/u/scratch/PCP4/PCP4kllt/project/
#SBATCH --partition medium
#SBATCH -t 30
#SBATCH -N 128
#SBATCH -n 8192
#SBATCH --overcommit
#SBATCH -o /gpfs/u/barn/PCP4/PCP4kllt/project/proj_128_08192.log

srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_GG.out --init 3 voronoi.00000.dat
srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_GG.out voronoi.00000.dat 10000 100

