CC=gcc
LIBS=-lm 
FLAGS= -O3 -pthread
qsort: qsort.c 
	$(CC) $(FLAGS) -o qsort qsort.c time.c $(LIBS)

clean:
	rm -f *.o *~ qsort time

