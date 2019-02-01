#!/bin/bash

a=$1

re='s/base/'$a'/g'
#re='s/1zdd/'$a'/g'
#re='s/max_iters: 5000/max_iters: 50000/'
#ree='s/max_evals: 500000/max_evals: 5000000/'
#re2='s/sade/sade_5k/'

for myfile in base*.yaml; do
    target=$(echo $myfile|sed -e $re)
    cp -v "$myfile" "$target"
    sed -i "$re" $target
done
