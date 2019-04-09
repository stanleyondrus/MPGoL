#!/bin/sh

#SBATCH --partition debug
#SBATCH --nodes 16
#SBATCH --time 60

srun -o gol_n16_r256_t1.log --ntasks 256 --overcommit assignment4-5.xl --threshold 0.25 --ticks 256 --write-uni --threads 1
srun -o gol_n16_r64_t4.log --ntasks 64 --overcommit assignment4-5.xl --threshold 0.25 --ticks 256 --write-uni --threads 4
srun -o gol_n16_r16_t16.log --ntasks 16 --overcommit assignment4-5.xl --threshold 0.25 --ticks 256 --write-uni --threads 16
srun -o gol_n16_r8_t32.log --ntasks 8 --overcommit assignment4-5.xl --threshold 0.25 --ticks 256 --write-uni --threads 32
srun -o gol_n16_r4_t64.log --ntasks 4 --overcommit assignment4-5.xl --threshold 0.25 --ticks 256 --write-uni --threads 64
