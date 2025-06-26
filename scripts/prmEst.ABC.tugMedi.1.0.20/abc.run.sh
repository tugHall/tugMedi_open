#!/bin/bash

cd $(dirname $0)
Rscript run.R >stdout.txt 2>stderr.txt
