/***************************************************************************/
/* Template for Asssignment 4/5 ********************************************/
/* Team Names Here              **(*****************************************/
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>

#include<clcg4.h>

#include<mpi.h>
#include<pthread.h>

// #define BGQ 1 // when running BG/Q, comment out when testing on mastiff

#ifdef BGQ
#include<hwi/include/bqc/A2_inlines.h>
#else
#define GetTimeBase MPI_Wtime
#endif

/***************************************************************************/
/* Defines *****************************************************************/
/***************************************************************************/

#define ALIVE 1
#define DEAD  0
#define NUM_TICKS 100
#define NUM_THREADS 3
typedef unsigned short int bool;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;

int uni_rows = 1024;
int uni_cols = 1024;
int gol_ticks = 10;

int mpi_myrank;
int mpi_ranks;

int rows_per_rank;

bool **universe;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

bool **alloc_bool_arr_2d(int rows, int cols);

void init_universe(bool **uni, int rows, int cols);

void sync_ghost_rows(bool **old_uni, int rows, int cols, int rank, int num_ranks);

void update_universe_state(bool **old_uni, int rows, int cols, int rank, int thread_id, double thresh);

bool life_lotto(Gen seed);

bool cell_next_state(int i, int j, bool **uni, int rows, int cols);

int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols);

void print_some_universe(bool **uni);

void free_bool_arr_2d(bool **uni, int rows);

void *threadFunction(void *arg);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]) {
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

    // Note, used the mpi_myrank to select which RNG stream to use.
    // You must replace mpi_myrank with the right row being used.
    // This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n",
           mpi_myrank, mpi_ranks, GenVal(mpi_myrank));

    MPI_Barrier(MPI_COMM_WORLD);

    // Start here
    rows_per_rank = uni_rows / mpi_ranks;

    // Initialize universe to all cells alive
    universe = alloc_bool_arr_2d(rows_per_rank + 2, uni_cols);
    init_universe(universe, rows_per_rank + 2, uni_cols);

    // TODO: Launch pthreads

    int thread_count = NUM_THREADS - 1;
    pthread_t tid[thread_count];
    int tid_index = 0;
    int rc;

    bool *arg = malloc(sizeof(bool));
    *arg = 69;

    int rows_per_thread = uni_rows / NUM_THREADS;

    for (int i = 0; i < thread_count; i++) {
        int *thread_id = malloc(sizeof(int));
        *thread_id = i;
        rc = pthread_create(&tid[tid_index], NULL, threadFunction, thread_id);
        if (rc != 0) { printf("Could not create thread"); }
        tid_index++;
    }

    // Run GoL simulation
    for (int t = 0; t < NUM_TICKS; t++) {
        sync_ghost_rows(universe, rows_per_rank + 2, uni_cols, mpi_myrank, mpi_ranks);
        if (t == 0) {
            update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, mpi_ranks, 0.70f);
        } else {
            update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, mpi_ranks, 0.0f);
        }
        if (mpi_myrank == 0 && t % 100 == 0) {
            printf("Tick: %d\n", t);
//            print_some_universe(universe);
        }
        // TODO: Write current state to file


//        int a = 4, b;
//        MPI_File fh;
//        MPI_File_open( MPI_COMM_WORLD, "myfile", MPI_MODE_RDWR, MPI_INFO_NULL, &fh ) ;
//        MPI_File_set_view( fh, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL ) ;
//        MPI_File_write_at(fh, 10, &a, 1, MPI_INT, &status ) ;
//        MPI_File_read_at(fh,  10, &b, 1, MPI_INT, &status ) ;
//        MPI_File_close( &cFile );

    }


    int *ret_val;
    for (int i = 0; i < tid_index; i++) {
        rc = pthread_join(tid[i], (void **) &ret_val);
        if (rc != 0) { printf("Could not join thread"); }
        printf("THREAD %ld: Thread [%ld] joined (returned %d)\n", pthread_self(), tid[i], *ret_val);
    }
    free(ret_val);


    // Clean up dynamically allocated variables
    free_bool_arr_2d(universe, rows_per_rank + 2);

    // END -Perform a barrier and then leave MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

void *threadFunction(void *arg) {
    int *thread_id = (int *) arg;

    printf("RANK %d THREAD %ld THREAD_ID %d\n", mpi_myrank, pthread_self(), *thread_id);

    for (int t = 0; t < NUM_TICKS; t++) {
        update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, *thread_id, (t == 0) ? 0.70f : 0.0f);
    }

//    free(arg);
    pthread_exit(thread_id);
}

void print_some_universe(bool **uni) {
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 50; j++) {
            if (uni[i][j] == ALIVE)
                printf("*");
            else
                printf(" ");
        }
        printf("\n");
    }
    printf("\n");
}

/*
 * Input:
 *  thresh - Value between 0.0 and 1.0. Changes the occurrence of random chance
 *          events where a cell is brought to live or killed.
 * We assume the universe is padded with 2 ghost rows, at the first and last indices.
 */
void update_universe_state(bool **old_uni, int rows, int cols, int rank, int thread_id, double thresh) {
    printf("UPDATING UNIVERSE AT RANK %d THREAD_ID %d THRESH %f\n", rank, thread_id, thresh);

    bool **new_uni = alloc_bool_arr_2d(rows, cols);
    for (int i = 0; i < rows; i++) {
        int global_row_ind = i + (rank * (rows - 2));
        for (int j = 0; j < cols; j++) {
            double random_chance = GenVal(global_row_ind);
            if (random_chance < thresh)
                new_uni[i][j] = life_lotto(global_row_ind);
            else
                new_uni[i][j] = cell_next_state(i, j, old_uni, rows, cols);
        }
    }

    // Copy new universe state into existing array
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            old_uni[i][j] = new_uni[i][j];
        }
    }

    free_bool_arr_2d(new_uni, rows);
}

bool life_lotto(Gen seed) {
    double chance = GenVal(seed);
    if (chance > 0.5f)
        return ALIVE;
    return DEAD;
}

void sync_ghost_rows(bool **old_uni, int rows, int cols, int rank, int num_ranks) {
    bool *tg_row = old_uni[0];
    bool *bg_row = old_uni[rows - 1];
    MPI_Request tg_req;
    MPI_Request bg_req;
    MPI_Status status1;
    MPI_Status status2;

    if (rank > 0) {
        // Recieve row from rank above, store in first "ghost" row
        MPI_Irecv(tg_row, cols, MPI_UNSIGNED_SHORT, rank - 1, 1, MPI_COMM_WORLD, &tg_req);

        // Send "first" row to rank above
        bool *f_row = old_uni[1];
        MPI_Request req1;
        MPI_Isend(f_row, cols, MPI_UNSIGNED_SHORT, rank - 1, 0, MPI_COMM_WORLD, &req1);
    }
    if (rank < num_ranks - 1) {
        // Recieve row from rank below, store in last "ghost" row
        MPI_Irecv(bg_row, cols, MPI_UNSIGNED_SHORT, rank + 1, 0, MPI_COMM_WORLD, &bg_req);

        // Send "last" row to rank below
        bool *l_row = old_uni[rows - 2];
        MPI_Request req2;
        MPI_Isend(l_row, cols, MPI_UNSIGNED_SHORT, rank + 1, 1, MPI_COMM_WORLD, &req2);
    }

    if (rank > 0) {
        MPI_Wait(&tg_req, &status1);
        // printf("Rank %d recieved top ghost row\n", rank);
    }
    if (rank < num_ranks - 1) {
        MPI_Wait(&bg_req, &status2);
        // printf("Rank %d recieved bottom ghost row\n", rank);
    }
}

int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols) {
    int count = 0;
    for (int i = i_in - 1; i < i_in + 2; i++) {
        for (int j = j_in - 1; j < j_in + 2; j++) {
            if (i > 0
                && j > 0
                && i < rows
                && j < cols
                && !(j == j_in && i == i_in)
                && uni[i][j] == ALIVE) {
                count++;
            }
        }
    }
    return count;
}

// 
bool cell_next_state(int i, int j, bool **uni, int rows, int cols) {
    // • Any live cell with fewer than two live neighbors dies, as if caused by under-population.
    // • Any live cell with more than three live neighbors dies, as if by over-population.
    // • Any live cell with two or three live neighbors lives on to the next generation.
    // • Any dead cell with exactly three live neighbors becomes a live cell, as if by reproduction
    int neighbors = neighbor_count(i, j, uni, rows, cols);
    int is_alive = uni[i][j];
    if (is_alive == ALIVE) {
        if (neighbors < 2 || neighbors > 3)
            return DEAD;
    } else {
        if (neighbors == 3)
            return ALIVE;
    }
    return is_alive;
}

void init_universe(bool **uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uni[i][j] = ALIVE;
        }
    }
}

bool **alloc_bool_arr_2d(int rows, int cols) {
    bool **uni = (bool **) calloc(rows, sizeof(bool *));
    for (int i = 0; i < rows; i++) {
        uni[i] = (bool *) calloc(cols, sizeof(bool));
    }
    return uni;
}

void free_bool_arr_2d(bool **uni, int rows) {
    for (int i = 0; i < rows; i++) {
        free(uni[i]);
    }
    free(uni);
}