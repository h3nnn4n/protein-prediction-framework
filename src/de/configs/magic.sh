#!/bin/bash

a=$1

re='s/base/'$a'/g'

for myfile in base*.yaml; do
    target=$(echo $myfile|sed -e $re)
    cp -v "$myfile" "$target"
    sed -i $re $a.yaml
done
