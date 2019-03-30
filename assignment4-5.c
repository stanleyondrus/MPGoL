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
#define bool unsigned short int

/***************************************************************************/
/* Global Vars *************************************************************/
/***************************************************************************/

double g_time_in_secs = 0;
double g_processor_frequency = 1600000000.0; // processing speed for BG/Q
unsigned long long g_start_cycles=0;
unsigned long long g_end_cycles=0;

int uni_rows = 1024;
int uni_cols = 1024;
int gol_ticks = 10;


/***************************************************************************/
/* Function Decs ***********************************************************/
/***************************************************************************/

void update_state(bool **rowin, bool **rowout);

void alloc_bool_arr_2d(bool ***uni, int rows, int cols);
void init_universe(bool ***uni, int rows, int cols);
void free_bool_arr_2d(bool ***uni, int rows);


/***************************************************************************/
/* Function: Main **********************************************************/
/***************************************************************************/

int main(int argc, char *argv[])
{
    //    int i = 0;
    int mpi_myrank;
    int mpi_ranks;
    // Example MPI startup and using CLCG4 RNG
    MPI_Init( &argc, &argv);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_ranks);
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_myrank);
    
    // Init 32,768 RNG streams - each rank has an independent stream
    InitDefault();
    
    // Note, used the mpi_myrank to select which RNG stream to use.
    // You must replace mpi_myrank with the right row being used.
    // This just show you how to call the RNG.    
    printf("Rank %d of %d has been started and a first Random Value of %lf\n", 
	   mpi_myrank, mpi_ranks, GenVal(mpi_myrank));
    
    MPI_Barrier( MPI_COMM_WORLD );

    // Start herevoid free_bool_arr_2d(bool ***uni, int rows)
    int rows_per_rank = uni_rows / mpi_ranks;

    // Allocate the universe array, split across MPI ranks (including "ghost" rows)
    bool **universe;
    alloc_bool_arr_2d(&universe, rows_per_rank + 2, uni_cols);

    // Set all 
    init_universe(&universe, rows_per_rank + 2, uni_cols);

    for (int t = 0; t < gol_ticks; t++) {
        universe_tick()
    }
    
    // Clean up dynamically allocated variables
    free_bool_arr_2d(&universe, rows_per_rank + 2);

    // END -Perform a barrier and then leave MPI
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize();
    return 0;
}

/***************************************************************************/
/* Other Functions - You write as part of the assignment********************/
/***************************************************************************/

void update_state(bool **rowin, bool **rowout) {

}

void init_universe(bool ***uni, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            (*uni)[i][j] = ALIVE;
        }
    }
}

void alloc_bool_arr_2d(bool ***uni, int rows, int cols) {
    *uni = (bool **)calloc(rows, sizeof(bool*));
    for (int i = 0; i < rows; i++) {
        (*uni)[i] = (bool *)calloc(cols, sizeof(bool));
    }
}

void free_bool_arr_2d(bool ***uni, int rows) {
    for (int i = 0; i < rows; i++) {
        free((*uni)[i]);
    }
    free(*uni);
}
