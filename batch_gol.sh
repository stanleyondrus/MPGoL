#!/bin/sh

#SBATCH --overcommit
#SBATCH --partition medium
#SBATCH --time 60
#SBATCH --nodes 128


# 4 Node
srun -o gol_n4_r64_t1.log --ntasks-per-node 64 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 1
srun -o gol_n4_r16_t4.log --ntasks-per-node 16 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 4
srun -o gol_n4_r4_t16.log --ntasks-per-node 4 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 16
srun -o gol_n4_r2_t32.log --ntasks-per-node 2 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 32
srun -o gol_n4_r1_t64.log --ntasks-per-node 1 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 64

# 16 Node
srun -o gol_n16_r64_t1.log --ntasks-per-node 64 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 1
srun -o gol_n16_r16_t4.log --ntasks-per-node 16 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 4
srun -o gol_n16_r4_t16.log --ntasks-per-node 4 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 16
srun -o gol_n16_r2_t32.log --ntasks-per-node 2 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 32
srun -o gol_n16_r1_t64.log --ntasks-per-node 1 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 64

# 64 Node
srun -o gol_n64_r64_t1.log --ntasks-per-node 64 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 1
srun -o gol_n64_r16_t4.log --ntasks-per-node 16 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 4
srun -o gol_n64_r4_t16.log --ntasks-per-node 4 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 16
srun -o gol_n64_r2_t32.log --ntasks-per-node 2 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 32
srun -o gol_n64_r1_t64.log --ntasks-per-node 1 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 64

# 128 Node
srun -o gol_n128_r64_t1.log --ntasks-per-node 64 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 1
srun -o gol_n128_r16_t4.log --ntasks-per-node 16 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 4
srun -o gol_n128_r4_t16.log --ntasks-per-node 4 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 16
srun -o gol_n128_r2_t32.log --ntasks-per-node 2 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 32
srun -o gol_n128_r1_t64.log --ntasks-per-node 1 assignment4-5.xl --threshold 0.25 --ticks 256 --write-heatmap --threads 64
