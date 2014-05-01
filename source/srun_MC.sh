#!/bin/bash
#
# srun command for CSCI-6360 group project.
# USAGE: /full/path/to/./q_MC.out [--help]
#                                 [--init dimension [outfile]]
#                                 [--nonstop dimension outfile steps [increment]]
#                                 [infile [outfile] steps [increment]]
#
#SBATCH --mail-type=END
#SBATCH --mail-user=kellet@rpi.edu
#SBATCH -D /gpfs/u/scratch/PCP4/PCP4kllt/project/
#SBATCH --partition small
#SBATCH -t 90
#SBATCH -N 64
#SBATCH -n 2048
#SBATCH --overcommit
#SBATCH -o /gpfs/u/barn/PCP4/PCP4kllt/project/proj_64_2048.log

# srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out --init 3 voronoi.0000.dat
# srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out voronoi.0000.dat 1500 50
srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out --nonstop 3 voronoi.0000.dat 1500 50
