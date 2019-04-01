# MPGoL
### Massively Parallel Game of Life with MPI, Pthreads, and Parallel I/O

CSCI 4320 - Parallel Computing and Programming

# Set Up
__This project requires a system which is capable of compiling and running MPI-based parallel programs.__

We included `sync_run_mastiff.sh`, which:
- Syncs local files to the mastiff server using _unison_
- Compiles the program using _make_
- Runs the program with _mpirun_
- Syncs remote files to local directory (in case your program writes to file)

You can tweak the number of MPI ranks to allocate by modifying the parameters of the `mpirun` command in the script.

# Objectives
### GoL with MPI
We achieved parallelization of the Game of Life:
1. Calculate rows_per_rank = rows_in_universe / ranks
2. Initialize universe, set all cells alive
3. Simulate the game, for each universe tick we update the state of each cell based on a set of rules which determines if the cell is alive or dead. We can also introduce randomn alive/dead events by tweaking the _random_chance_ threshold.
### GoL with MPI and Pthreads
To be implemented...
### Parallel File I/O
To be implemented...
### Heatmap
To be implemented...

# Experiments
### Strong Scaling with Parallel I/O
__Run tests on these configurations__

| nodes | 0 threads/rank 64 MPI ranks/node | 4 threads/rank 16 MPI ranks/node | 16 threads/rank 4 MPI ranks/node | 32 threads/rank 2 MPI ranks/node | 64 threads/rank 1 MPI rank/node |
|-------|----------------------------------|----------------------------------|----------------------------------|----------------------------------|---------------------------------|
| 4     |                                  |                                  |                                  |                                  |                                 |
| 16    |                                  |                                  |                                  |                                  |                                 |
| 64    |                                  |                                  |                                  |                                  |                                 |
| 128   |                                  |                                  |                                  |                                  |                                 |
### Parallel I/O and Heatmap of Final Universe State
