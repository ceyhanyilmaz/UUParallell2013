
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <math.h>

#define DEBUG 0

void quickSort( int[], int, int, int);
int partition( int[], int, int);

void quickSort(int a[], int l, int r, int level) {
   int j;

   if( l < r ) {
      j = partition( a, l, r);

      if (level != 0) {
        quickSort( a, l, j-1, level-1);
        quickSort( a, j+1, r, level-1);
      } else {
        #pragma omp task
        quickSort( a, l, j-1, level-1);
        #pragma omp task
        quickSort( a, j+1, r, level-1);
      }
   }
}

int partition(int a[], int l, int r) {
  int pivot, i, j, t;
  pivot = a[l];
  i = l; j = r+1;

  while( 1) {
    do ++i; while( a[i] <= pivot && i <= r );
    do --j; while( a[j] > pivot );

    if( i >= j ) break;

    t = a[i]; 
    a[i] = a[j]; 
    a[j] = t;
  }
  
  t = a[l]; 
  a[l] = a[j]; 
  a[j] = t;

  return j;
}

void main (int argc, char *argv[]) {
  int elements_maximum, thread_maximum,t,j,i;
  int *elements;
  int levels = 6;
  if (argc != 3) {
    printf("Example of usage: ./qsort <elements> <threads>\n");
    return;
  }

  elements_maximum = atoi(argv[1]);
  thread_maximum = atoi(argv[2]);

  omp_set_num_threads(thread_maximum);
  elements = (int *)malloc(elements_maximum*sizeof(int));  

  for(j = 0; j <= elements_maximum; j++) {
    elements[j] = (rand() % 999 +1);
  }

  if (DEBUG) {
    printf("\n\nUnsorted array is:  ");
    for(i = 0; i < elements_maximum; ++i) {
      printf(" %d ", elements[i]);
    }
    printf("\n");
  }

  t = timer();

  #pragma omp parallel
  {
    #pragma omp single nowait
    {
        quickSort( elements, 0, elements_maximum-1, levels);
    }
  }

  t = timer() - t;
  printf("Completed %fs\n", t/1000000.0);

  if (DEBUG) {
    printf("\n\nSorted array is:  ");
    for(i = 0; i < elements_maximum; ++i) { 
      printf(" %d ", elements[i]);
    }
  } 

}
