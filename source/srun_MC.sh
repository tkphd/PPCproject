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
#SBATCH -t 30
#SBATCH -N 64
#SBATCH -n 2048
#SBATCH --overcommit
#SBATCH -o /gpfs/u/barn/PCP4/PCP4kllt/project/prMC_64_2048.log

# srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out --init 3 voronmc.000.dat
# srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out voronmc.000.dat 200 200
srun --runjob-opts="--mapping TEDCBA" /gpfs/u/barn/PCP4/PCP4kllt/project/source/./q_MC.out --nonstop 3 voronmc.000.dat 200 200 2
