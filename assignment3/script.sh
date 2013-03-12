#!/bin/bash

file='ass3_gram.output'

> file

echo "gram started"
for n in 1000 2000 3000 4000 5000
do
	echo "n = $n started"
	echo "n = $n" >> $file
	for t in 1 2 4 8 16 32 64
	do
		echo "t = $t" >> $file

		./gram $n $t >> $file
		./gram $n $t >> $file
		./gram $n $t >> $file
		./gram $n $t >> $file
		echo "t = $t completed"
	done
done
echo "gram completed"
