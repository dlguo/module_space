#!/bin/sh

SESS=$1
for i in $(ls $DATA/results_SIFT2/); do
  qsub -q standby -v SUBJ=${i},SESS=${SESS} get_nets.sub
done
