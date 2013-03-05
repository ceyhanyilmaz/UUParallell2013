#!/bin/bash

file='ass1.output'

> file

echo "fox started"
for n in 100 400 600 800 1000 1200
do
	echo "n = $n started"
	echo "n = $n" >> $file
	for t in 1 4 9 16
	do
		echo "t = $t" >> $file
        mpirun -hostfile nodes -mca plm_rsh_agent rsh -np $t ass1 $n ass1 >> $file
        mpirun -hostfile nodes -mca plm_rsh_agent rsh -np $t ass1 $n ass1 >> $file
        mpirun -hostfile nodes -mca plm_rsh_agent rsh -np $t ass1 $n ass1 >> $file
        mpirun -hostfile nodes -mca plm_rsh_agent rsh -np $t ass1 $n ass1 >> $file
		echo "t = $t completed"
	done
done
echo "fox completed"
