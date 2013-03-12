
CC=gcc
LIBS=-lm 
CFLAGS    =  -O3 -fopenmp

gram: gram.c 
	$(CC) $(CFLAGS) -o gram gram.c time.c $(LIBS)

qsort: qsort.c 
	$(CC) $(CFLAGS) -o qsort qsort.c time.c $(LIBS)

clean:
	rm -f *.o *~ gram time

