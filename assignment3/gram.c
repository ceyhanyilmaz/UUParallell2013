#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>
 
double **V,**Q;

double vecNorm(double *,int );
double scalarProd(double *,double *,int);


int main(int argc, char *argv[]) {
   int n,i,j,k,time,threads;
   double t,sigma,temp_norm,norm;

   if (argc != 3) {
      printf("Example of usage: ./gram <elements> <threads>\n");
      return;
   }

   n = atoi(argv[1]);
   threads = atoi(argv[2]);

   omp_set_num_threads(threads);
  
   //Allocate and fill vectors
   V = (double **)malloc(n*sizeof(double *));
   Q = (double **)malloc(n*sizeof(double *));

   for(i=0; i<n; i++) {
      V[i] = (double *)malloc(n*sizeof(double));
      Q[i] = (double *)malloc(n*sizeof(double));
   }

   for (i = 0; i<n; i++) {
      for(j=0; j<n; j++) {
        V[i][j] = rand() % 5 + 1;
      }
   }
    
   time=timer();

   for(i=0; i<n; i++) {
      temp_norm = vecNorm(V[i],n);
      for (k=0; k<n; k++) {
         Q[i][k] = V[i][k]/temp_norm;
      }

      #pragma omp parallel for private (k,sigma)
      for(j=i+1; j<n; j++) {
         sigma = scalarProd(Q[i],V[j],n);
         for(k=0;k<n;k++) {
            V[j][k] -=sigma*Q[i][k];
         }
     }
}
   time=timer()-time;
   printf("Elapsed time: %f \n",time/1000000.0);
   
   return 0;
}

double vecNorm(double *vec,int n) {
   int i;
   double local_norm = 0;
   
   for(i=0;i<n;i++) {
      local_norm+= (vec[i]*vec[i]);
   }

   return sqrt(local_norm);
}

double scalarProd(double *a,double *b,int n) {
   int i;
   double scalar =0.0;
   
   for(i=0;i<n;i++) {
      scalar +=a[i]*b[i];
   }
   
   return scalar;
}
