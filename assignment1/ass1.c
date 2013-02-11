
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define N_MAX 10000000

int main(int argc, char *argv[]) {
  int rank, size,n,row,col,count,blocklen,stride,seed,myid,nproc,i,j,subrank,data,mydata;
  double *A, *B, *tmp;
  MPI_Status status;
  MPI_Datatype newtype;
  MPI_Comm proc_grid, proc_col, proc_row;
  int coords[2],pos[2],reorder=1,ndim=2,dims[2]={0,0},periods[2]={0,0};

  seed = time(NULL);
  srand(seed);
  
  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  /* TODO: CHECK THAT size IS SQUARE */
  n=size; /* Size of matrices depend on number of processes (?) */
  //nx=atoi(argv[1]); ny=atoi(argv[2]); 					/* Set size of matrices */
	A=(double *)calloc(n*n,sizeof(double)); 		/* Matrix A */
	B=(double *)calloc(n*n,sizeof(double)); 		/* Matrix B */	
  tmp=(double *)calloc(n*n,sizeof(double));		/* temporary matrix */
	
	if (rank==0){
    printf("Matrix A on proc 0\n");
    for (row=0; row<n;row++){
      for (col=0; col<n;col++){
        A[row*n+col]= (rand() % 9 +1);
        B[row*n+col]= (rand() % 9 +1);
        printf("%d ", (int)A[row*n+col]);
      }
      printf("\n");
    }
    printf("\n");
  } 

	/* Create a virtual 2D-grid topology */
  MPI_Dims_create(nproc,ndim,dims);
  MPI_Cart_create(MPI_COMM_WORLD,ndim,dims,periods,reorder,&proc_grid);
  MPI_Comm_rank(proc_grid,&myid);    /* Note: use proc_grid */

	MPI_Cart_coords(proc_grid,myid,ndim,coords); 
	   
	 /* Create a communicator for each row */
  MPI_Comm_split(proc_grid,coords[0],coords[1],&proc_row);
  MPI_Comm_rank(proc_row,&subrank);
  
   /* Create a communicator for each column */
  MPI_Comm_split(proc_grid,coords[1],coords[0],&proc_col);
  MPI_Comm_rank(proc_col,&subrank);



	/* TODO: Distribute matrices A and B so that each process gets A_local and B_local */
	


	/* Do the FOX */
	
	/* Collect submatrices C_local (Allgather?) */
	
	



  count=n/2; blocklen=n/2; stride=n;
  MPI_Type_vector(count,blocklen,stride,MPI_DOUBLE,&newtype);
  MPI_Type_commit(&newtype);

  /* Send last column, notice the message length = 1 ! */
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

  free(A);
  free(B);
  free(tmp);
  MPI_Type_free(&newtype);
  MPI_Finalize(); 

  return 0;
}
