/***************************************************************************/
/* Asssignment 4/5  ********************************************************/
/* Judy Fan         ********************************************************/
/* Stanley Ondrus   ********************************************************/
/* Sean Rice        ********************************************************/

// Compile with:
// BGQ: make b
// Mastiff: make all

// Run with:
// sbatch --partition debug --nodes 4 --time 60 ./run_n4.sh
// sbatch --partition debug --nodes 16 --time 60 ./run_n16.sh
// sbatch --partition small --nodes 64 --time 60 ./run_n64.sh
// sbatch --partition medium --nodes 128 --time 60 ./run_n128.sh
/***************************************************************************/

/***************************************************************************/
/* Includes ****************************************************************/
/***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>
#include<math.h>
#include<unistd.h>

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
#define true 1
#define false 0
typedef unsigned short int bool;

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles = 0;
unsigned long long g_end_cycles = 0;

int num_threads = 16;
int uni_rows = 4096;
int uni_cols = 4096;
int heatmap_rows = 0;
int heatmap_cols = 0;
float threshold = 0.50f;
int gol_ticks = 1;
bool write_heatmap = false;
bool write_universe = false;

int mpi_myrank;
int mpi_ranks;

int uni_rows_per_rank;
int heatmap_rows_per_rank;
int heatmap_row_inc;
int heatmap_col_inc;
int rows_per_thread;

bool **universe;
bool **new_uni;
int **heatmap;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t presync_barrier, barrier;

struct prog_params {
    int uni_rows;
    int uni_cols;
    int heatmap_rows;
    int heatmap_cols;
    int ticks;
    int threads;
    float threshold;
    bool write_universe;
    bool write_heatmap;
};

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

// Parameter parsing
void print_params(struct prog_params params);
struct prog_params parse_args(int argc, char **argv);

// Universe ops
void init_universe(bool **uni, int rows, int cols);
void sync_ghost_rows(bool **old_uni, int rows, int cols, int rank, int num_ranks);
void update_universe_state(bool **old_uni, int rows, int cols, int rank, int thread_id, double thresh);
bool cell_next_state(int i, int j, bool **uni, int rows, int cols);
bool life_lotto(Gen seed);
int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols);
void print_some_universe(bool **uni, int rows, int cols);
void write_universe_to_file(char *filename, bool **uni, int rows, int cols, int rank);

// Heatmaps ops
void calculate_heatmap();
void print_some_heatmap(int **uni, int rows, int cols);
void write_heatmap_to_file(char *filename, int **uni, int rows, int cols, int rank);
int sum_matrix(bool **uni, int x_start, int x_stop, int y_start, int y_stop);

// Alloc/free 2d array funcs
bool **alloc_bool_arr_2d(int rows, int cols);
int **alloc_int_arr_2d(int rows, int cols);
void free_bool_arr_2d(bool **uni, int rows);
void free_int_arr_2d(int **uni, int rows);

// Pthread funcs
void *threadFunction(void *arg);

// Helper funcs
int mod(int a, int b);
void strreverse(char *begin, char *end);
void itoa(int value, char *str, int base);

/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]) {
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    // Parse command-line args
    struct prog_params params = parse_args(argc, argv);
    num_threads = params.threads;
    uni_rows = params.uni_rows;
    uni_cols = params.uni_cols;
    heatmap_rows = params.heatmap_rows;
    heatmap_cols = params.heatmap_cols;
    gol_ticks = params.ticks;
    threshold = params.threshold;
    write_heatmap = params.write_heatmap;
    write_universe = params.write_universe;

    // Print rank and thread configuration
    if (mpi_myrank == 0) {
        printf("Ranks: %d\n", mpi_ranks);
        print_params(params);
    }

    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize global variables
    uni_rows_per_rank = uni_rows / mpi_ranks;
    heatmap_rows_per_rank = heatmap_rows / mpi_ranks;
    heatmap_col_inc = uni_cols / heatmap_cols;
    heatmap_row_inc = uni_rows / heatmap_rows;

    // Initialize universe to all cells alive
    universe = alloc_bool_arr_2d(uni_rows_per_rank + 2, uni_cols);
    new_uni = alloc_bool_arr_2d(uni_rows_per_rank + 2, uni_cols);
    init_universe(universe, uni_rows_per_rank + 2, uni_cols);

    heatmap = alloc_int_arr_2d(heatmap_rows_per_rank, heatmap_cols);

    int thread_ids[num_threads - 1];
    pthread_t tid[num_threads - 1];
    int rc;

    rows_per_thread = uni_rows_per_rank / num_threads;

    pthread_barrier_init(&presync_barrier, NULL, num_threads);

    for (int i = 0; i < num_threads - 1; i++) {
        thread_ids[i] = i + 1;
        rc = pthread_create(&tid[i], NULL, threadFunction, &thread_ids[i]);
        if (rc != 0) { printf("Could not create thread"); }
    }

    g_start_cycles = GetTimeBase();

    // Run GoL simulation
    for (int t = 0; t < gol_ticks; t++) {
        sync_ghost_rows(universe, uni_rows_per_rank + 2, uni_cols, mpi_myrank, mpi_ranks);
        // printf("Sync ghost rows\n");

        pthread_barrier_wait(&presync_barrier);
        // printf("post barrier");
//        pthread_barrier_init(&barrier, NULL, thread_count + 1);
//        pthread_barrier_wait(&barrier);

        // ranks are considered thread 0
        update_universe_state(universe, uni_rows_per_rank + 2, uni_cols, mpi_myrank, 0, (t == 0) ? threshold : 0.0f);


        pthread_barrier_wait(&presync_barrier);
        if (write_heatmap == true) {
            calculate_heatmap();
        }


        if (mpi_myrank == 0) {
            printf("Tick: %d\n", t);
//            print_some_universe(universe, uni_rows_per_rank + 2, uni_cols);
//            print_universe_to_file("file.txt", universe, 32, 32);
//            print_some_heatmap(heatmap, uni_rows_per_rank / HEATMAP_WINDOW_SIZE, HEATMAP_SIZE);
        }
    }

    g_end_cycles = GetTimeBase();
    if (mpi_myrank == 0) {
        g_time_in_secs = ((double) (g_end_cycles - g_start_cycles)) / g_processor_frequency;
        printf("Main execution time: %lf\n", g_time_in_secs);
    }

    for (int i = 0; i < num_threads - 1; i++) {
        int *ret_val;
        rc = pthread_join(tid[i], (void **) &ret_val);
        if (rc != 0) { printf("Could not join thread"); }
//        printf("THREAD %ld: Thread [%ld] joined (returned %d)\n", pthread_self(), tid[i], *ret_val);
    }
    // free(ret_val);

    MPI_Barrier(MPI_COMM_WORLD);

    if (write_heatmap == true) {
        if (mpi_myrank == 0) {
            printf("Writing heatmap to file.\n");
        }
        write_heatmap_to_file("heatmap_out.ibin", heatmap, heatmap_rows_per_rank, heatmap_cols, mpi_myrank);
    }
    if (write_universe == true) {
        if (mpi_myrank == 0) {
            printf("Writing universe to file.\n");
        }
        g_start_cycles = GetTimeBase();
        write_universe_to_file("uni_out.usbin", universe, uni_rows_per_rank + 2, uni_cols, mpi_myrank);
        g_end_cycles = GetTimeBase();
        if (mpi_myrank == 0) {
            g_time_in_secs = ((double) (g_end_cycles - g_start_cycles)) / g_processor_frequency;
            printf("Parallel I/O execution time: %lf\n", g_time_in_secs);
        }
    }
    // Clean up dynamically allocated variables
    free_bool_arr_2d(universe, uni_rows_per_rank + 2);
    free_int_arr_2d(heatmap, heatmap_rows_per_rank);

    // END - Perform a barrier and then leave MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Argument Parsing ********************************************************/
/***************************************************************************/

struct prog_params parse_args(int argc, char **argv) {
    struct prog_params params;
    params.threads = 1;
    params.uni_rows = 32768;
    params.uni_cols = 32768;
    params.heatmap_cols = 1024;
    params.heatmap_rows = 1024;
    params.ticks = 256;
    params.threshold = 0.50f;
    params.write_heatmap = false;
    params.write_universe = false;

    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp("--threads", argv[i]) == 0) {
                params.threads = atoi(argv[i + 1]);
            }
            if (strcmp("--ticks", argv[i]) == 0) {
                params.ticks = atoi(argv[i + 1]);
            }
            if (strcmp("--uni-size", argv[i]) == 0) {
                params.uni_rows = atoi(argv[i + 1]);
                params.uni_cols = atoi(argv[i + 2]);
            }
            if (strcmp("--heatmap-size", argv[i]) == 0) {
                params.heatmap_rows = atoi(argv[i + 1]);
                params.heatmap_cols = atoi(argv[i + 2]);
            }
            if (strcmp("--write-uni", argv[i]) == 0) {
                params.write_universe = true;
            }
            if (strcmp("--write-heatmap", argv[i]) == 0) {
                params.write_heatmap = true;
            }
            if (strcmp("--threshold", argv[i]) == 0) {
                params.threshold = atof(argv[i + 1]);
            }
        }
    }

    return params;
}

void print_params(struct prog_params params) {
    printf("Threads per Rank: %d\n", params.threads);
    printf("Universe Size: %d by %d\n", params.uni_rows, params.uni_cols);
    printf("Heatmap Size: %d by %d\n", params.heatmap_rows, params.heatmap_cols);
    printf("Ticks: %d\n", params.ticks);
    printf("Threshold: %f\n", params.threshold);
}

/***************************************************************************/
/* Pthread *****************************************************************/
/***************************************************************************/

void *threadFunction(void *arg) {
    int *thread_id = (int *) arg;

//    printf("RANK %d THREAD %ld THREAD_ID %d\n", mpi_myrank, pthread_self(), *thread_id);

//    pthread_mutex_lock(&mutex);
//    pthread_mutex_unlock(&mutex);


    for (int t = 0; t < gol_ticks; t++) {
        // Wait for the MPI_ranks to sync
        pthread_barrier_wait(&presync_barrier);
        // printf("RANK %d THREAD %ld THREAD_ID %d iteration: %d\n", mpi_myrank, pthread_self(), *thread_id, t);
        update_universe_state(universe, uni_rows_per_rank + 2, uni_cols, mpi_myrank, *thread_id, (t == 0) ? threshold : 0.0f);
        // Wait for all threads to finish syncing before printing
        pthread_barrier_wait(&presync_barrier);
    }

//    free(arg);
    pthread_exit(thread_id);
}

/***************************************************************************/
/* Universe Operations *****************************************************/
/***************************************************************************/

void sync_ghost_rows(bool **old_uni, int rows, int cols, int rank, int num_ranks) {
    bool *tg_row = old_uni[0];
    bool *bg_row = old_uni[rows - 1];
    MPI_Request tg_req;
    MPI_Request bg_req;
    MPI_Status status1;
    MPI_Status status2;

    // Recieve row from rank above, store in first "ghost" row
    MPI_Irecv(tg_row, cols, MPI_UNSIGNED_SHORT, mod(rank - 1, num_ranks), 1, MPI_COMM_WORLD, &tg_req);
    // Recieve row from rank below, store in last "ghost" row
    MPI_Irecv(bg_row, cols, MPI_UNSIGNED_SHORT, mod(rank + 1, num_ranks), 0, MPI_COMM_WORLD, &bg_req);

    // Send "first" row to rank above
    bool *f_row = old_uni[1];
    MPI_Request req1;
    MPI_Isend(f_row, cols, MPI_UNSIGNED_SHORT, mod(rank - 1, num_ranks), 0, MPI_COMM_WORLD, &req1);

    // Send "last" row to rank below
    bool *l_row = old_uni[rows - 2];
    MPI_Request req2;
    MPI_Isend(l_row, cols, MPI_UNSIGNED_SHORT, mod(rank + 1, num_ranks), 1, MPI_COMM_WORLD, &req2);


    MPI_Wait(&tg_req, &status1);
    // printf("Rank %d recieved top ghost row\n", rank);
    MPI_Wait(&bg_req, &status2);
    // printf("Rank %d recieved bottom ghost row\n", rank);
}

void update_universe_state(bool **old_uni, int rows, int cols, int rank, int thread_id, double thresh) {
    /*
    * Input:
    *  thresh - Value between 0.0 and 1.0. Changes the occurrence of random chance
    *          events where a cell is brought to live or killed.
    * We assume the universe is padded with 2 ghost rows, at the first and last indices.
    */

    // Allocate a new universe to temporarily store next universe state
    // Row bounds for a thread to update
    int start = (thread_id * rows_per_thread) + 1;
    int end = start + rows_per_thread - 1;
    if (thread_id == num_threads - 1) {
        end = rows - 1;
    }

    // Update all cells
    for (int i = start; i < end && i < rows; i++) {
        int global_row_ind = i + (rank * (rows - 2));
        for (int j = 0; j < cols; j++) {
            double random_chance = GenVal(global_row_ind);
            if (random_chance <= thresh) {
                new_uni[i][j] = life_lotto(global_row_ind);
            } else
                new_uni[i][j] = cell_next_state(i, j, old_uni, rows, cols);
        }
    }

    // Copy new universe state into existing array
    for (int i = start; i < end && i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            old_uni[i][j] = new_uni[i][j];
        }
    }
}

void init_universe(bool **uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            uni[i][j] = ALIVE;
        }
    }
}

// Coin flip: return ALIVE/DEAD with 50% chance each
bool life_lotto(Gen seed) {
    double chance = GenVal(seed);
    if (chance > 0.5f)
        return ALIVE;
    return DEAD;
}


bool cell_next_state(int i, int j, bool **uni, int rows, int cols) {
    /* Rules:
    *
    * - Any live cell with fewer than two live neighbors dies, as if caused
    * by under-population.
    * 
    * - Any live cell with more than three live neighbors dies, as if by
    * over-population.
    * 
    * - Any live cell with two or three live neighbors lives on to the next
    * generation.
    * 
    * - Any dead cell with exactly three live neighbors becomes a live cell, as
    * if by reproduction
    */
    int neighbors = neighbor_count(i, j, uni, rows, cols);
    int is_alive = uni[i][j];
    if (is_alive == ALIVE) {
        if (neighbors < 2 || neighbors > 3) {
            return DEAD;
        }
    } else {
        if (neighbors == 3) {
            return ALIVE;
        }
    }
    return is_alive;
}

int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols) {
    int count = 0;
    // Check all eight neighbors of input cell
    for (int i = i_in - 1; i < i_in + 2; i++) {
        for (int j = j_in - 1; j < j_in + 2; j++) {
            // Verify index is:
            //  - inside of array
            //  - index is not "center": (i_in, j_in)
            if (i > 0 && j > 0 && i < rows && j < cols &&
                !(j == j_in && i == i_in) && uni[i][j] == ALIVE) {
                count++;
            }
        }
    }
    return count;
}

void print_some_universe(bool **uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (uni[i][j] == ALIVE) {
                printf("*");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }
    printf("\n");
}


void write_universe_to_file(char *filename, bool **uni, int rows, int cols, int rank) {
    MPI_Status status;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    for (int i = 1; i < rows - 1; i++) {
        // printf("%d\n", rows-2);
        int global_row_ind = (i - 1) + (rank * (rows-2));
        int offset = (global_row_ind * cols * 2);
        // printf("Rank %d writing row %d with offset %d.\n", rank, i, offset);
        MPI_File_write_at(fh, offset, uni[i], cols, MPI_UNSIGNED_SHORT, &status);
    }

    MPI_File_close(&fh);
}


/***************************************************************************/
/* Heatmap Operations ******************************************************/
/***************************************************************************/

void calculate_heatmap() {
    for (int i = 0; i < heatmap_rows_per_rank; i++) {
        for (int j = 0; j < heatmap_cols; j++) {
            int start_x = heatmap_row_inc * i;
            int start_y = heatmap_col_inc * j;
            int stop_x = start_x + heatmap_row_inc;
            int stop_y = start_y + heatmap_col_inc;
            heatmap[i][j] += sum_matrix(universe, start_x, stop_x, start_y, stop_y);
        }
    }
}

int sum_matrix(bool **uni, int x_start, int x_stop, int y_start, int y_stop) {
    int ret_val = 0;
    for (int i = x_start; i > 0 && i < uni_rows_per_rank + 1 && i < x_stop; ++i) {
        for (int j = y_start; j < y_stop; ++j) {
            ret_val += uni[i][j];
        }
    }
    return ret_val;
}

void print_some_heatmap(int **uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %d ", uni[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void write_heatmap_to_file(char *filename, int **uni, int rows, int cols, int rank) {
    MPI_Status status;
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);

    for (int i = 0; i < rows; i++) {
        int global_row_ind = i + (rank * (rows));
        int offset = (global_row_ind * cols * 4);
        // printf("Rank %d writing row %d with offset %d.\n", rank, i, offset);
        MPI_File_write_at(fh, offset, uni[i], cols, MPI_INT, &status);
    }

    MPI_File_close(&fh);
}


/***************************************************************************/
/* Alloc/Free 2D Array Funcs ***********************************************/
/***************************************************************************/

int **alloc_int_arr_2d(int rows, int cols) {
    int **uni = (int **) calloc(rows, sizeof(int *));
    for (int i = 0; i < rows; i++) {
        uni[i] = (int *) calloc(cols, sizeof(int));
    }
    return uni;
}

void free_int_arr_2d(int **uni, int rows) {
    for (int i = 0; i < rows; i++) {
        free(uni[i]);
    }
    free(uni);
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


/***************************************************************************/
/* Helper Functions ********************************************************/
/***************************************************************************/

int mod(int a, int b) {
    // We can get the correct modulus for negative values using this function
    int r = a % b;
    return r < 0 ? r + b : r;
}

void strreverse(char *begin, char *end) {
    char aux;
    while (end > begin) {
        aux = *end, *end-- = *begin, *begin++ = aux;
    }
}

void itoa(int value, char *str, int base) {
    static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char *wstr = str;
    int sign;

    // Validate base
    if (base < 2 || base > 35) {
        *wstr = '\0';
        return;
    }

    // Take care of sign
    if ((sign = value) < 0) value = -value;

    // Conversion. Number is reversed.
    do *wstr++ = num[value % base]; while (value /= base);
    if (sign < 0) *wstr++ = '-';
    *wstr = '\0';

    // Reverse string
    strreverse(str, wstr - 1);
}