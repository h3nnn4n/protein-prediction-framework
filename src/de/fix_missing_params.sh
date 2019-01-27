#!/bin/bash

for file in stats__*.dat
do
  base=$(basename $file .dat)
  target=$(echo $base | sed -e 's/stats/parameters/').yaml
  touch $target
done
