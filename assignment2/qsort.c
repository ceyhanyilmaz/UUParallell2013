// quickSort.c
#include <stdio.h>
#include <stdlib.h> 
#include <pthread.h>
#include <time.h>
#include <math.h>

#define DEBUG 0

void qsort_serial( int l, int r);
void *qsort_parallel(void *arg);
int partition_serial( int l, int r);

pthread_mutex_t mutexsort;

struct thread_data {
  int start_index;
  int end_index;
};

int thread_count = 1;
int thread_maximum = 0;
int elements_maximum = 0;
int *elements;

void *qsort_parallel(void *arg) {
  struct thread_data index_left, index_right;
  int start_index, end_index, pivot;

  void *status;
  pthread_t thread1, thread2;
  
  struct thread_data *args;
  args = (struct thread_data *) arg;
  start_index = args->start_index;
  end_index = args->end_index;

  if (start_index < end_index) {
    
    if (thread_count < thread_maximum-2) {
      pivot = partition_serial(start_index, end_index);

      index_left.start_index = start_index;
      index_left.end_index = pivot-1;

      index_right.start_index = pivot+1;
      index_right.end_index = end_index;

      pthread_mutex_lock(&mutexsort);
      thread_count += 2;
      pthread_mutex_unlock(&mutexsort);

      pthread_create(&thread1, NULL, qsort_parallel,(void *) &index_left);
      pthread_create(&thread2, NULL, qsort_parallel,(void *) &index_right);

      pthread_join(thread1, &status);
      pthread_join(thread2, &status);

    } else {
      qsort_serial(start_index, end_index);
    }
  }

  pthread_exit(NULL);
}

void qsort_serial( int l, int r) {
   int pivot;

   if (l < r)  {
      pivot = partition_serial(l, r);

      qsort_serial(l, pivot-1);
      qsort_serial(pivot+1, r);
   }
}

int partition_serial( int l, int r) {
  int pivot, i, j, t; 
  pivot = elements[l];
  i = l; 
  j = r+1;
    
  while (1) {
    do 
      ++i; 
    while (elements[i] <= pivot && i <= r);

    do 
      --j; 
    while (elements[j] > pivot);

    if( i >= j ) break;

    t = elements[i]; 
    elements[i] = elements[j]; 
    elements[j] = t;
  }
  
  t = elements[l]; 
  elements[l] = elements[j]; 
  elements[j] = t;
  
  return j;
}

void main (int argc, char *argv[]) {
  int seed, rank, nthreads, t, i, j, rc;
  struct thread_data index_left;
  void *status;
  int pivot;
  int ttime;

  if (argc != 3) {
    printf("Example of usage: ./qsort <elements> <threads>\n");
    return;
  }

  elements_maximum = atoi(argv[1]);
  thread_maximum = atoi(argv[2]);

  ttime = timer();
  seed = time(NULL);
  srand(seed);

  pthread_t thread;
  pthread_attr_t attr;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&mutexsort, NULL);

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
  
  index_left.start_index = 0;
  index_left.end_index =  elements_maximum-1;

  rc = pthread_create(&thread, &attr, qsort_parallel, (void *) &index_left);

  pthread_join(thread, &status);

  if (DEBUG) {
    printf("\n\nSorted array is:  ");
    for(i = 0; i < elements_maximum; ++i) {
      printf(" %d ", elements[i]);
    }
    printf("\n");
  }

  ttime = timer()-ttime;
  printf("\nTime: %f, Threads used: %d\n\n", ttime/1000000.0, thread_count);

  pthread_mutex_destroy(&mutexsort);
  pthread_attr_destroy(&attr);
  free(elements);
  
  pthread_exit(NULL);
}

