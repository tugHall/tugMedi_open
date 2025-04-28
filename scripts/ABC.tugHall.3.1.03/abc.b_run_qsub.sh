#!/bin/bash

RUN_LIST=./run_list.txt

count=1
while read script; do
  qsub -N "simJob_$count" "$script"
  count=$(($count + 1))
done < ${RUN_LIST}
