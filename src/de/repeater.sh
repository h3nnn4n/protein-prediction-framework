#!/bin/bash

ns=$1
file=$2

for i in `seq $ns`
do
    python3.5 main.py $file
done
