CC=gcc
LIBS=-lm 
FLAGS= -O3 -pthread
qsort: qsort.c 
	$(CC) $(FLAGS) -o qsort qsort.c $(LIBS)

qsort2: qsort2.c 
	$(CC) $(FLAGS) -o qsort2 qsort2.c time.c $(LIBS)

clean:
	rm -f *.o *~ qsort qsort2 time

