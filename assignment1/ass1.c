#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define N_MAX 10000000

typedef struct {
    int       size;           /* # of processes   		*/
    MPI_Comm  proc_grid;      /* Grid communicator 		*/
    MPI_Comm  proc_row;       /* Row communicator  		*/
    MPI_Comm  proc_col;  			/* Column communicator  */
    int       length;    			/* Length of grid       */
    int       rowrank;    		/* My row number        */
    int       colrank;    		/* My column number     */
    int       gridrank;   		/* Grid rank     				*/
} GRID_INFO_T;

void setup_grid(GRID_INFO_T *grid){
	int rank;
	int coords[2],pos[2],reorder=1,ndim=2,dims[2],periodic[2]={1,1};
		
	MPI_Comm_size(MPI_COMM_WORLD, &(grid->size)); 	 /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 					 /* Get my number                */
	
	grid->length = (int)sqrt((double)grid->size);
  assert(grid->length*grid->length == grid->size); /* Check to see if size is a good square number */
  
  dimensions[2] = {grid->length,grid->length};
  
  MPI_Dims_create(grid->size,ndim,dims);
  MPI_Cart_create(MPI_COMM_WORLD,ndim,dims,periodic,reorder,&proc_grid); /* Create grid */
  MPI_Comm_rank(proc_grid,&gridrank);    					 /* Distribute grid ranks */

	MPI_Cart_coords(proc_grid,gridrank,ndim,coords); /* Gives coordinates for gridrank */
	   
	 /* Create a communicator for each row */
  MPI_Comm_split(proc_grid,coords[0],coords[1],&proc_row);
  MPI_Comm_rank(proc_row,&rowrank);
  
   /* Create a communicator for each column */
  MPI_Comm_split(proc_grid,coords[1],coords[0],&proc_col);
  MPI_Comm_rank(proc_col,&colrank);
}

void fill_matrix(double* matrix){
	//printf("Matrix:\n");
    for (row=0; row<n;row++){
      for (col=0; col<n;col++){
        matrix[row*n+col]= (rand() % 9 +1);
        //printf("%d ", (int)matrix[row*n+col]);
      }
     // printf("\n");
    }
   // printf("\n");
}


int main(int argc, char *argv[]) {
	GRID_INFO_T grid;
  double *local_A, *local_B, *local_C;
  double *A, *B, *C;
  int n,seed,count,blocklen,stride;
  MPI_Datatype submatrix;
  
  seed = time(NULL);
  srand(seed);
	n=argv[1]; 	


  
/* Fill and print matrix A and B */
	
	if (rank==0){
		A=(double *)calloc(n*n,sizeof(double)); 	/* Matrix A */
	  B=(double *)calloc(n*n,sizeof(double)); 	/* Matrix B */	
    C=(double *)calloc(n*n,sizeof(double));		/* Matrix C */
		fill_matrix(A);
		fill_matrix(B);
  } 
  
	MPI_Init(&argc, &argv); 		/* Initialize MPI               */
	setup_grid(&grid);
	
	count=blocklen=n/grid->length; /* Number of elements in matrix divided by length of grid(sqrt of processes) */
	stride=n;											 /* Number of elements in one row/column -> n */
 
  MPI_Type_vector(count,blocklen,stride,MPI_DOUBLE,&submatrix);
  MPI_Type_commit(&submatrix);

	if (rank == 0) {
    MPI_Send(&A[n/2*n+n], 1, newtype, 1, 111, MPI_COMM_WORLD);

  } else if (rank==1) {
    MPI_Recv(&A[n/2*n], 1, newtype, 0, 111, MPI_COMM_WORLD, &status);
    printf("Matrix A on proc 1\n");
    for (row=0; row<n; row++)
    {
       for (col=0; col<n; col++)
          printf("%4d ", (int)A[row*n+col]);
       printf("\n"); 
    }
  }
	/* TODO: Distribute matrices A and B so that each process gets A_local and B_local */
	


	/* Do the FOX */
	
	/* Collect submatrices C_local (Allgather?) */
	
	



  count=n/2; blocklen=; stride=n;
  MPI_Type_vector(count,blocklen,stride,MPI_DOUBLE,&newtype);
  MPI_Type_commit(&newtype);

  /* Send last column, notice the message length = 1 ! */
  

  free(A);
  free(B);
  free(C);
  MPI_Type_free(&newtype);
  MPI_Finalize(); 

  return 0;
  /*
  int rank, size,n,row,col,count,blocklen,stride,seed,myid,nproc,i,j,colrank,rowrank,data,mydata, elemsBlock;

  MPI_Status status;
  MPI_Datatype newtype;
  MPI_Comm proc_grid, proc_col, proc_row;
  int coords[2],pos[2],reorder=1,ndim=2,dims[2],periods[2]={0,0}; */
}




