/*
 * Assignment: Fox implementation 
 * Author: Anders Hassis, Jill Karlsson & Andreas Moregård
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#define N_MAX 10000000
#define DEBUG 0

typedef struct GRID_INFO_T {
    int       size;             /* # of processes           */
    MPI_Comm  proc_grid;        /* Grid communicator        */
    MPI_Comm  proc_row;         /* Row communicator         */
    MPI_Comm  proc_col;         /* Column communicator  */
    int       length;           /* Length of grid       */
    int       rowrank;          /* My row number        */
    int       colrank;          /* My column number     */
    int       gridrank;         /* Grid rank                    */
} GRID_INFO_T;

void printMatrix(double *m, int dims) {
    int x,y;
    for (x = 0; x < dims; x++) {
        for (y = 0; y < dims; y++) {
            printf("%3.0f ", m[x*dims+y]);
        }
        printf("\n");
    }
    printf("\n");
}

void PrintGridInfo (GRID_INFO_T *grid) {
    printf ("--------------------\n");
    printf ("Number of Processes is %d\n", grid->size);
    printf ("Grid Comm Identifier is %d\n", grid->proc_grid);
    printf ("Row Comm Identifier is %d\n", grid->proc_row);
    printf ("Column Comm Identifier is %d\n", grid->proc_col);
    printf ("Grid Order is %d\n", grid->length);
    printf ("Current Process Coordinates are (%d, %d)\n", grid->rowrank, grid->colrank);
    printf ("Process rank in Grid is %d\n", grid->gridrank); 
    printf ("--------------------\n");
}

void setup_grid(GRID_INFO_T *grid) {
    int rank, coords[2], pos[2], reorder = 1, ndim = 2, dims[2], periodic[2] = {0,0};
        
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);            /* Get my number                */
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->size));    /* Get the number of processors */
    
    grid->length = (int)sqrt((double)grid->size);

    assert(grid->length*grid->length == grid->size); /* Check to see if size is a good square number */
    dims[0] = dims[1] = grid->length;

    MPI_Dims_create(grid->size, ndim, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periodic, reorder, &(grid->proc_grid));  /* Create grid */
    MPI_Comm_rank(grid->proc_grid, &(grid->gridrank));                                   /* Distribute grid ranks */

    MPI_Cart_coords(grid->proc_grid,grid->gridrank,ndim,coords);                         /* Gives coordinates for gridrank */

     /* Create a communicator for each row */
    MPI_Comm_split(grid->proc_grid,coords[0],coords[1],&(grid->proc_row));
    MPI_Comm_rank(grid->proc_row,&(grid->rowrank));

    /* Create a communicator for each column */
    MPI_Comm_split(grid->proc_grid,coords[1],coords[0],&(grid->proc_col));
    MPI_Comm_rank(grid->proc_col,&(grid->colrank));
}

void fill_matrix(double* matrix,int n){
    if (DEBUG) { printf("Matrix:\n"); }
    int row,col;
    for (row=0; row<n;row++) {
      for (col=0; col<n;col++) {
        matrix[row*n+col]= (rand() % 2 +1);
        if (DEBUG) { printf("%d ", (int)matrix[row*n+col]); }
      }
     if (DEBUG) { printf("\n"); }
    }
   if (DEBUG) { printf("\n"); }
}

// Multiply two matrices (m1 and m2) and put result in m3
void multiplyLocal(int n, double *m1, double *m2, double *m3) {
    int i,j,k;
     for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                m3[(i*n)+j] += m1[(i*n)+k]*m2[(k*n)+j];
            }
        }
    }
}
 
int Fox(int elems, GRID_INFO_T *grid, double *A, double *B, double *C) {
    double *tempA, *tempB; 
    MPI_Status status, status2;
    MPI_Request request,request2;

    tempA = (double *)calloc(elems*elems,sizeof(double));  
    tempB = (double *)calloc(elems*elems,sizeof(double));  
    int stage, root;

    for (stage = 0; stage<grid->length; stage++) {
        root = (grid->colrank+stage)%(grid->length);

        if (root == grid->rowrank) {
            for (int i = 0; i < elems*elems; i++) {
                tempA[i] = A[i];
            }
        }
            
        MPI_Bcast(tempA, elems*elems, MPI_DOUBLE, root, grid->proc_row);

        // Send B
        if (grid->colrank == 0) {
            MPI_Isend(B, elems*elems, MPI_DOUBLE, (grid->length)-1, 113, grid->proc_col, &request);
        } else {
            MPI_Isend(B, elems*elems, MPI_DOUBLE, (grid->colrank)-1, 113, grid->proc_col, &request);
        }

        if (grid->colrank == (grid->length)-1) {
            MPI_Irecv(tempB, elems*elems, MPI_DOUBLE, 0, 113, grid->proc_col, &request);
        } else {
            MPI_Irecv(tempB, elems*elems, MPI_DOUBLE, (grid->colrank)+1, 113, grid->proc_col, &request);
        }

        multiplyLocal(elems, tempA, B, C);

        MPI_Wait(&request, &status);
        for (int i=0; i < elems*elems; i++) {
            B[i] = tempB[i];
        }
    }
}

void collectMatrix(double *global_C, double *local_C, GRID_INFO_T *grid, int n, MPI_Datatype *subMatrix) {
    MPI_Status status;
    MPI_Request request;
    int elems = n/grid->length;
    int comm_rank, position;

    MPI_Isend(local_C, elems*elems, MPI_DOUBLE, 0, 0, grid->proc_grid, &request);

    if (grid->gridrank == 0) {
        for(int i=0; i < grid->size; i++) {
            int col = floor(i/grid->length)*elems;
            int row = (i % grid->length)*elems;
            int position = row*n+col;

            MPI_Irecv(&global_C[position], 1, *subMatrix, i, 0, grid->proc_grid, &request);
            MPI_Wait(&request, &status);
        }
    }
}

int main(int argc, char *argv[]) {
    double startTime = 0.0, endTime = 0.0, finish = 0.0;
    double *global_A, *global_B, *global_C;
    double *A, *B, *C;
    int n, seed, count, blocklen, stride;

    GRID_INFO_T grid;
    MPI_Status status;
    MPI_Datatype submatrix;
  
    seed = time(NULL);
    srand(seed);
    n = atoi(argv[1]); 
  
    MPI_Init(&argc, &argv);             /* Initialize MPI */
    MPI_Barrier(MPI_COMM_WORLD);

    setup_grid(&grid);

    stride = n;                         /* Number of elements in one row/column -> n */
    count = blocklen = n/grid.length;   /* Number of elements in matrix divided by length of grid(sqrt of processes) */
    A = (double *)calloc(count*count,sizeof(double));
    B = (double *)calloc(count*count,sizeof(double));
    C = (double *)calloc(count*count,sizeof(double));
    
    /* Fill and print matrix A and B */
    if (grid.gridrank == 0) {
        global_A = (double *)calloc(n*n,sizeof(double));    /* Matrix A */
        global_B = (double *)calloc(n*n,sizeof(double));    /* Matrix B */  
        global_C = (double *)calloc(n*n,sizeof(double));    /* Matrix C */

        fill_matrix(global_A, n);
        fill_matrix(global_B, n);
    }

    MPI_Type_vector(count, blocklen, stride, MPI_DOUBLE, &submatrix);
    MPI_Type_commit(&submatrix);

    startTime = MPI_Wtime();

    int elems = n/grid.length;
    /* Distribute matrices A and B so that each process gets A_local and B_local */
    if (grid.gridrank == 0) {
        int i;
        for(i=1; i < grid.size; i++) {
            int col = floor(i/grid.length)*elems;
            int row = (i % grid.length)*elems;

            MPI_Send(&global_A[row*n+col], 1, submatrix, i, 111, grid.proc_grid);
            MPI_Send(&global_B[row*n+col], 1, submatrix, i, 112, grid.proc_grid);
        }
        int j,k;
        for (j = 0; j < n/grid.length; j++) {
            for (k = 0; k < n/grid.length; k++) {
                A[j*(n/grid.length)+k] = global_A[j*n+k];
                B[j*(n/grid.length)+k] = global_B[j*n+k];
            }
        }
    } else {  
        MPI_Recv(A, elems*elems, MPI_DOUBLE, 0, 111, grid.proc_grid, &status);
        MPI_Recv(B, elems*elems, MPI_DOUBLE, 0, 112, grid.proc_grid, &status);
    }

    /* Do the FOX */
    Fox(elems, &grid, A, B, C);

    /* Syncronize the calculations */
    MPI_Barrier(MPI_COMM_WORLD);
    collectMatrix(global_C, C, &grid, n, &submatrix);

    if (grid.gridrank == 0) {
        endTime = MPI_Wtime();
        finish = endTime-startTime;
        printf("Execution time: %.10f\n", finish);
    }

    /* Print our finalized C matrix */ 
    if (grid.gridrank == 0 && DEBUG) {
        printf("Finished product:\n");
        int frow,fcol;

        for (frow=0; frow<n;frow++) {
            for (fcol=0; fcol<n;fcol++) {
                printf("%3.0f ", global_C[frow*n+fcol]);
            }
            printf("\n");
        }
        printf("\n");   
    }


    free(A);
    free(B);
    free(C);

   
    MPI_Comm_free(&(grid.proc_col));
    MPI_Comm_free(&(grid.proc_row));
    MPI_Comm_free(&(grid.proc_grid));

    MPI_Type_free(&submatrix);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize(); 

  return 0;
}
