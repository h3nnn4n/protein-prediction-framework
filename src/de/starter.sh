#!/bin/bash

a=$1

for i in `seq 8`
do
	./runner.sh $i $a
done
