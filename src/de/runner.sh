#!/bin/bash

a=$1
b=$2


if [ $a -eq 1 ]; then
for i in `seq 10`; do
		x='configs/'$b'_lhs_10_mc_pop.yaml'
		python3.5 main.py $x
		x='configs/'$b'_lhs_10.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 2 ]; then
for i in `seq 10`; do
		x='configs/'$b'_lhs_20_mc_pop.yaml'
		python3.5 main.py $x
		x='configs/'$b'_lhs_20.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 3 ]; then
for i in `seq 10`; do
		x='configs/'$b'_lhs_5_mc_pop.yaml'
		python3.5 main.py $x
		x='configs/'$b'_lhs_5.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 4 ]; then
for i in `seq 10`; do
		x='configs/'$b'_ls_best_100.yaml'
		python3.5 main.py $x
		x='configs/'$b'_ls_best_250.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 5 ]; then
for i in `seq 10`; do
		x='configs/'$b'_ls_best_500.yaml'
		python3.5 main.py $x
		x='configs/'$b'_ls_pop_100.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 6 ]; then
for i in `seq 10`; do
		x='configs/'$b'_ls_pop_250.yaml'
		python3.5 main.py $x
		x='configs/'$b'_ls_pop_500.yaml'
		python3.5 main.py $x
		x='configs/'$b'_reset_1000.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 7 ]; then
for i in `seq 10`; do
		x='configs/'$b'_reset_250.yaml'
		python3.5 main.py $x
		x='configs/'$b'_reset_500.yaml'
		python3.5 main.py $x
		x='configs/'$b'_reset_750.yaml'
		python3.5 main.py $x
done
elif [ $a -eq 8 ]; then
for i in `seq 10`; do
		x='configs/'$b'_stage0.yaml'
		python3.5 main.py $x
		x='configs/'$b'_v.yaml'
		python3.5 main.py $x
done
fi
