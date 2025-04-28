#!/bin/bash
export LANG=

NJOB='90%'
LOG=log/simJob_$(date '+%Y%m%d_%H%M%S').txt
RUN_LIST=./run_list.txt

mkdir -p log
parallel \
    -j $NJOB \
    --joblog ${LOG} \
    --progress \
    -a ${RUN_LIST}
