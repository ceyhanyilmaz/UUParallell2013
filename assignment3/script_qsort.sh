#!/bin/bash

file='ass3_qsort.output'

> file

echo "qsort started"
for n in 100000 1000000 10000000
do
	echo "n = $n started"
	echo "n = $n" >> $file
	for t in 1 2 4 8 16 32 64
	do
		echo "t = $t" >> $file

		./qsort $n $t >> $file
		./qsort $n $t >> $file
		./qsort $n $t >> $file
		./qsort $n $t >> $file
		echo "t = $t completed"
	done
done
echo "qsort completed"
