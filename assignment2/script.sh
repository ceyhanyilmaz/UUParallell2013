#!/bin/bash

file='ass2.output'

> file

echo "qsort started"
for n in 1000 10000 100000 1000000 10000000
do
	echo "n = $n started"
	echo "n = $n" >> $file
	for t in 1 5 10 25 50 75 100 150
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
