#!/bin/bash

sed -i 's|/rshare1/ZETTAI_path_WA_slash_home_KARA||g' run_list.txt
sed -i 's|/abc.run.sh||g' run_list.txt

RUN_LIST=./run_list.txt
shell=$SHELL

count=1
while read script; do
  cd $script
  qsub \
  -N tugHall_$count \
  -S $shell \
  -cwd \
  -o $script \
  -e $script \
  $script/abc.run_qsub.sh
  count=$(($count + 1))
done < ${RUN_LIST}
