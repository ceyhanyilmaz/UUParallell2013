#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#define N_MAX 10000000

typedef struct GRID_INFO_T {
    int       size;             /* # of processes   		*/
    MPI_Comm  proc_grid;        /* Grid communicator 		*/
    MPI_Comm  proc_row;         /* Row communicator  		*/
    MPI_Comm  proc_col;  		/* Column communicator  */
    int       length;    		/* Length of grid       */
    int       rowrank;    		/* My row number        */
    int       colrank;    		/* My column number     */
    int       gridrank;   		/* Grid rank     				*/
} GRID_INFO_T;

void PrintGridInfo (GRID_INFO_T *grid) {
    printf ("--------------------\n");
	printf ("Number of Processes is %d\n", grid->size);
	printf ("Grid Comm Identifier is %d\n", grid->proc_grid);
	printf ("Row Comm Identifier is %d\n", grid->proc_row);
	printf ("Column Comm Identifier is %d\n", grid->proc_col);
	printf ("Grid Order is %d\n", grid->length);
	printf ("Current Process Coordinates are (%d, %d)\n",
            grid->rowrank, grid->colrank);
	printf ("Process rank in Grid is %d\n", grid->gridrank); 
    printf ("--------------------\n");
}

void setup_grid(GRID_INFO_T *grid){
	int rank;
	int coords[2],pos[2],reorder=1,ndim=2,dims[2],periodic[2]={1,1};
		
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 			 /* Get my number                */
	MPI_Comm_size(MPI_COMM_WORLD, &(grid->size)); 	 /* Get the number of processors */
	

	grid->length = (int)sqrt((double)grid->size);
    assert(grid->length*grid->length == grid->size); /* Check to see if size is a good square number */
    dims[0] = dims[1] = grid->length;

    MPI_Dims_create(grid->size, ndim, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periodic, reorder, &(grid->proc_grid)); /* Create grid */
    MPI_Comm_rank(grid->proc_grid, &(grid->gridrank));    					 /* Distribute grid ranks */

    MPI_Cart_coords(grid->proc_grid,grid->gridrank,ndim,coords); /* Gives coordinates for gridrank */

	 /* Create a communicator for each row */
    MPI_Comm_split(grid->proc_grid,coords[0],coords[1],&(grid->proc_row));
    MPI_Comm_rank(grid->proc_row,&(grid->rowrank));

    /* Create a communicator for each column */
    MPI_Comm_split(grid->proc_grid,coords[1],coords[0],&(grid->proc_col));
    MPI_Comm_rank(grid->proc_col,&(grid->colrank));

}

void fill_matrix(double* matrix,int n){
	printf("Matrix:\n");
    int row,col;
    //int n = (int)sqrt((double)sizeof(matrix));

    for (row=0; row<n;row++){
      for (col=0; col<n;col++){
        matrix[row*n+col]= (rand() % 9 +1);
        printf("%d ", (int)matrix[row*n+col]);
      }
     printf("\n");
    }
   printf("\n");
}


void multiplyLocal(int n, double *m1, double *m2, double *m3) {
    int i,j,k;
     for (i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                m3[(j*n)+i] += m1[(k*n)+i]*m2[(j*n)+k];
            }

            printf("C_(%d,%d) = %f\n", i,j,m3[(j*n)+i]);
        }
    }
}
 
int Fox(int n, MPI_Datatype *submatrix, GRID_INFO_T *grid, double *A, double *B, double *C) {
    double *tempMatrix; 
    MPI_Status status;
    
    tempMatrix = (double *)calloc(n*n,sizeof(double));	
    int i, root;
    for (i = 0; i<grid->length; i++) {
        root = i; //(grid->rowrank+i)%(grid->length)
        printf("Root is %d\n", root);
        
        if (root == grid->colrank) {
            MPI_Bcast(A, n/grid->length, MPI_DOUBLE, root, grid->proc_row);
            multiplyLocal(n/grid->length, A, B, C);
        } else {
            MPI_Bcast(tempMatrix, n/grid->length, MPI_DOUBLE, root, grid->proc_row);
            multiplyLocal(n/grid->length, tempMatrix, B, C);
        }

        MPI_Sendrecv_replace(B, 1, MPI_DOUBLE, (grid->rowrank+grid->length+1)%grid->length, 0, 
            (grid->rowrank+grid->length-1)%grid->length, 0, grid->proc_row, &status);
    }
}


int main(int argc, char *argv[]) {
	GRID_INFO_T grid;
    double *local_A, *local_B, *local_C;
    double *A, *B, *C;
    int n, seed, count, blocklen, stride;
        MPI_Status status;
    MPI_Datatype submatrix;
  
    seed = time(NULL);
    srand(seed);
	n = atoi(argv[1]); 
  
	MPI_Init(&argc, &argv); 		 /* Initialize MPI */
	setup_grid(&grid);

	stride=n;					     /* Number of elements in one row/column -> n */
	count=blocklen=n/grid.length;    /* Number of elements in matrix divided by length of grid(sqrt of processes) */

    /* Fill and print matrix A and B */
	if (grid.gridrank == 0) {
        A = (double *)calloc(n*n,sizeof(double)); 	/* Matrix A */
	    B = (double *)calloc(n*n,sizeof(double)); 	/* Matrix B */	
        C = (double *)calloc(n*n,sizeof(double));	/* Matrix C */

		fill_matrix(A, n);
		fill_matrix(B, n);

    } else {
        A = (double *)calloc(count,sizeof(double));
        B = (double *)calloc(count,sizeof(double));
        C = (double *)calloc(n*n,sizeof(double));
    }
    //PrintGridInfo(&grid);

    MPI_Type_vector(count,blocklen,stride,MPI_DOUBLE,&submatrix);
    MPI_Type_commit(&submatrix);

    int elems = n/grid.length;

	/* Distribute matrices A and B so that each process gets A_local and B_local */
	if (grid.gridrank == 0) {
        int i;

        for(i=1; i<grid.size; i++) {
            int col = floor(i/grid.length)*elems;
            int row = (i % grid.length)*elems;

            MPI_Send(&A[row*n+col], 1, submatrix, i, 111, grid.proc_grid);
            MPI_Send(&B[row*n+col], 1, submatrix, i, 111, grid.proc_grid);
        }
    } else {  
        MPI_Recv(&A[0], 1, submatrix, 0, 111, grid.proc_grid, &status);
        MPI_Recv(&B[0], 1, submatrix, 0, 111, grid.proc_grid, &status);

        int i, j;
        int row = grid.gridrank % grid.length;


        for (i=0; i<elems; i++) {
            for (j=0; j<elems; j++) {
                   if (grid.gridrank == 1) {
                printf("Element: %4d on processor %d ", (int)A[i*n+j],grid.gridrank);
}
            }
            printf("\n"); 
    }
  }

	/* Do the FOX */
    Fox(n, &submatrix, &grid, A, B, C);

	/* Collect submatrices C_local (Allgather?) */
	MPI_Barrier(MPI_COMM_WORLD);
    
    int i;
    for (i = 0

  /* Send last column, notice the message length = 1 ! */
  if (grid.gridrank == 0) {
    free(A);
    free(B);
    free(C);
  }

  MPI_Type_free(&submatrix);

  MPI_Finalize(); 

  return 0;
  /*
  int rank, size,n,row,col,count,blocklen,stride,seed,myid,nproc,i,j,colrank,rowrank,data,mydata, elemsBlock;

  MPI_Status status;
  MPI_Datatype newtype;
  MPI_Comm proc_grid, proc_col, proc_row;
  int coords[2],pos[2],reorder=1,ndim=2,dims[2],periods[2]={0,0}; */
}




