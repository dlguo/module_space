#!/bin/sh

FUNC=$1
SESS=$2
WINSIZE=$3
CUTOFF=$4
if [ ${FUNC} = "dtw" ]; then
  qsub -q standby -v FUNC=${FUNC},SUBJ="NA",SESS=${SESS},WINSIZE=${WINSIZE},CUTOFF=${CUTOFF} main.sub
else
  for i in $(ls $DATA/results_SIFT2/); do
    qsub -q standby -v FUNC=${FUNC},SUBJ=${i},SESS=${SESS},WINSIZE=${WINSIZE},CUTOFF=${CUTOFF} main.sub
  done
fi