/***************************************************************************/
/* Asssignment 4/5  ********************************************************/
/* Judy Fan         ********************************************************/
/* Stanley Ondrus   ********************************************************/
/* Sean Rice        ********************************************************/
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
#define NUM_TICKS 256
#define NUM_THREADS 16
#define HEATMAP_SIZE 1024
#define HEATMAP_WINDOW_SIZE 32
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
int gol_ticks = 1;

int mpi_myrank;
int mpi_ranks;

int rows_per_rank; // TODO: probably take this one out
int rows_per_thread;

bool **universe;
int **heatmap;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t presync_barrier, barrier;

/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

bool **alloc_bool_arr_2d(int rows, int cols);

int **alloc_int_arr_2d(int rows, int cols);

void init_universe(bool **uni, int rows, int cols);
void sync_ghost_rows(bool **old_uni, int rows, int cols, int rank, int num_ranks);
void update_universe_state(bool **old_uni, int rows, int cols, int rank, int thread_id, double thresh);
bool life_lotto(Gen seed);
bool cell_next_state(int i, int j, bool **uni, int rows, int cols);
int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols);
void print_some_universe(bool **uni, int rows, int cols);

void print_some_heatmap(int **uni, int rows, int cols);

void free_bool_arr_2d(bool **uni, int rows);
void *threadFunction(void *arg);
int mod(int a, int b);
void write_universe_to_file(const char *filename, int **uni, int rows, int cols, int rank);

void print_universe_to_file(const char *filename, bool **uni, int rows, int cols);

int sum_matrix(bool **uni, int x_start, int x_stop, int y_start, int y_stop);

void calculate_heatmap();


/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[]) {
    // Example MPI startup and using CLCG4 RNG
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);

    // Parse command-line args
    if (argc >= 2) {
        num_threads = atoi(argv[1]);
    }

    // Print rank and thread configuration
    if (mpi_myrank == 0)
        printf("Ranks: %d\nThreads per Rank: %d\n", mpi_ranks, num_threads);

    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();

    // Note, used the mpi_myrank to select which RNG stream to use.
    // You must replace mpi_myrank with the right row being used.
    // This just show you how to call the RNG.    
    // printf("Rank %d of %d has been started and a first Random Value of %lf\n",
    //        mpi_myrank, mpi_ranks, GenVal(mpi_myrank));

    MPI_Barrier(MPI_COMM_WORLD);

    // Start here
    rows_per_rank = uni_rows / mpi_ranks;

    // Initialize universe to all cells alive
    universe = alloc_bool_arr_2d(rows_per_rank + 2, uni_cols);
    init_universe(universe, rows_per_rank + 2, uni_cols);
    printf("Rank %d created local universe with %d rows.\n", mpi_myrank, rows_per_rank + 2);

    heatmap = alloc_int_arr_2d(rows_per_rank / HEATMAP_WINDOW_SIZE, uni_cols / HEATMAP_WINDOW_SIZE);
    printf("Rank %d created local heatmap (%d, %d).\n", mpi_myrank, rows_per_rank / HEATMAP_WINDOW_SIZE, uni_cols / HEATMAP_WINDOW_SIZE);

    // TODO: Launch pthreads

    int thread_count = num_threads - 1;
    pthread_t tid[thread_count];
    int tid_index = 0;
    int rc;

    rows_per_thread = uni_rows / (num_threads * mpi_ranks);

    pthread_barrier_init(&presync_barrier, NULL, num_threads);

    for (int i = 0; i < thread_count; i++) {
        int *thread_id = malloc(sizeof(int));
        *thread_id = i + 1;
        rc = pthread_create(&tid[tid_index], NULL, threadFunction, thread_id);
        if (rc != 0) { printf("Could not create thread"); }
        tid_index++;
    }

    // Run GoL simulation
    for (int t = 0; t < NUM_TICKS; t++) {
        sync_ghost_rows(universe, rows_per_rank + 2, uni_cols, mpi_myrank, mpi_ranks);
        pthread_barrier_wait(&presync_barrier);
//        pthread_barrier_init(&barrier, NULL, thread_count + 1);
//        pthread_barrier_wait(&barrier);
        if (t == 0) {
            // ranks are considered thread 0
            update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, 0, 0.70f);
//            update_universe_state(universe, rows_per_thread + 2, uni_cols, mpi_myrank, 0, 0.70f);
        } else {
            update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, 0, 0.0f);
//            update_universe_state(universe, rows_per_thread + 2, uni_cols, mpi_myrank, 0, 0.0f);
        }
        pthread_barrier_wait(&presync_barrier);
        calculate_heatmap();

        if (mpi_myrank == 0 && t % 10 == 0) {
            printf("Tick: %d\n", t);
//            print_some_universe(universe, rows_per_rank + 2, uni_cols);
//            print_universe_to_file("file.txt", universe, 32, 32);
//            print_some_heatmap(heatmap, rows_per_rank / HEATMAP_WINDOW_SIZE, HEATMAP_SIZE);
        }

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
//        printf("THREAD %ld: Thread [%ld] joined (returned %d)\n", pthread_self(), tid[i], *ret_val);
    }
    free(ret_val);

    write_universe_to_file("/home/parallel/spring-2019/rices/hw4/uni_out.ibin", heatmap,rows_per_rank / HEATMAP_WINDOW_SIZE, uni_cols / HEATMAP_WINDOW_SIZE, mpi_myrank);
    // Clean up dynamically allocated variables
    free_bool_arr_2d(universe, rows_per_rank + 2);

    // END -Perform a barrier and then leave MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Function Implementations ************************************************/
/***************************************************************************/

// We can get the correct modulus for negative values using this function
int mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

void calculate_heatmap() {
    for (int i = 0; i < rows_per_rank / HEATMAP_WINDOW_SIZE; i++) {
        for (int j = 0; j < uni_cols / HEATMAP_WINDOW_SIZE; j++) {
            int start_x = HEATMAP_WINDOW_SIZE * i;
            int start_y = HEATMAP_WINDOW_SIZE * j;
            int stop_x = start_x + HEATMAP_WINDOW_SIZE;
            int stop_y = start_y + HEATMAP_WINDOW_SIZE;
            heatmap[i][j] += sum_matrix(universe, start_x, stop_x, start_y, stop_y);
        }
    }
}

int sum_matrix(bool **uni, int x_start, int x_stop, int y_start, int y_stop) {
    int ret_val = 0;
    for (int i = x_start; i < x_stop; ++i) {
        for (int j = y_start; j < y_stop; ++j) {
            ret_val += uni[i][j];
        }
    }
    return ret_val;
}

void *threadFunction(void *arg) {
    int *thread_id = (int *) arg;

//    printf("RANK %d THREAD %ld THREAD_ID %d\n", mpi_myrank, pthread_self(), *thread_id);

//    pthread_mutex_lock(&mutex);
//    pthread_mutex_unlock(&mutex);


    for (int t = 0; t < NUM_TICKS; t++) {
        // Wait for the MPI_ranks to sync
        pthread_barrier_wait(&presync_barrier);
//        printf("RANK %d THREAD %ld THREAD_ID %d iteration: %d\n", mpi_myrank, pthread_self(), *thread_id, t);
        update_universe_state(universe, rows_per_rank + 2, uni_cols, mpi_myrank, *thread_id, (t == 0) ? 0.70f : 0.0f);
        // Wait for all threads to finish syncing before printing
        pthread_barrier_wait(&presync_barrier);
    }

//    free(arg);
    pthread_exit(thread_id);
}

void strreverse(char* begin, char* end) {
	
	char aux;
	
	while(end>begin) {
		aux=*end, *end--=*begin, *begin++=aux;
    }
	
}

void itoa(int value, char* str, int base) {
	
	static char num[] = "0123456789abcdefghijklmnopqrstuvwxyz";
	
	char* wstr=str;
	
	int sign;
	
	// Validate base	
	if (base<2 || base>35){ *wstr='\0'; return; }
	
	// Take care of sign
	if ((sign=value) < 0) value = -value;

	// Conversion. Number is reversed.
	do *wstr++ = num[value%base]; while(value/=base);
	if(sign<0) *wstr++='-';
	*wstr='\0';
	
	// Reverse string
	strreverse(str,wstr-1);
}

void write_universe_to_file(const char *filename, int **uni, int rows, int cols, int rank) {
    MPI_Status status;
    MPI_File fh;
    MPI_File_open( MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh ) ;

    for (int i = 0; i < rows; i++) {
        int global_row_ind = i + (rank * (rows));
        int offset = (global_row_ind * cols * 4);
        printf("Rank %d writing row %d with offset %d.\n", rank, i, offset);
        MPI_File_write_at(fh, offset, uni[i], cols, MPI_INT, &status );
    }

    MPI_File_close( &fh );
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

void print_some_heatmap(int **uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %d ", uni[i][j]);
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
//    printf("UPDATING UNIVERSE AT RANK %d THREAD_ID %d THRESH %f\n", rank, thread_id, thresh);

    // Allocate a new universe to temporarily store next universe state
    bool **new_uni = alloc_bool_arr_2d(rows, cols);

    // Row bounds for a thread to update
    int start = thread_id * (rows_per_thread) + 1;
    int end = (start + ((rows - 2) / num_threads));
    if (thread_id == num_threads - 1) {
        end = rows - 1;
    }

    // Update all cells in rows
//    printf("For loop i=%d; i<%d\n", start, end);
    for (int i = start; i < end; i++) {
        int global_row_ind = i + (rank * (rows - 2));
//        printf("rank %d, rows %d, thread_id %d, global row_ind = %d \n", rank, rows, thread_id, global_row_ind);
//        printf("thread: %d, rows: %d, i: %d, global_row_ind: %d\n", thread_id, rows, i, global_row_ind);
        for (int j = 0; j < cols; j++) {
            double random_chance = GenVal(global_row_ind);
            if (random_chance < thresh) {
                new_uni[i][j] = life_lotto(global_row_ind);
//                printf("global row index: %d, random chance\n", global_row_ind);
            } else
                new_uni[i][j] = cell_next_state(i, j, old_uni, rows, cols);
        }
    }

    // Copy new universe state into existing array
    for (int i = start; i < end; i++) {
        for (int j = 0; j < cols; j++) {
            old_uni[i][j] = new_uni[i][j];
        }
    }

    free_bool_arr_2d(new_uni, rows);
}

// Coin flip: return ALIVE/DEAD with 50% chance each
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

int neighbor_count(int i_in, int j_in, bool **uni, int rows, int cols) {
    int count = 0;
    // Check all eight neighbors of input cell
    for (int i = i_in - 1; i < i_in + 2; i++) {
        for (int j = j_in - 1; j < j_in + 2; j++) {
            // Verify index is:
            //  - inside of array
            //  - index is not "center"
            if (i > 0 && j > 0 && i < rows && j < cols && 
                !(j == j_in && i == i_in) && uni[i][j] == ALIVE) {
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

int **alloc_int_arr_2d(int rows, int cols) {
    int **uni = (int **) calloc(rows, sizeof(int *));
    for (int i = 0; i < rows; i++) {
        uni[i] = (int *) calloc(cols, sizeof(int));
    }
    return uni;
}

void free_bool_arr_2d(bool **uni, int rows) {
    for (int i = 0; i < rows; i++) {
        free(uni[i]);
    }
    free(uni);
}
