#!/bin/sh

if [ "$1" == "gen" ]; then
  SESS=$2
  WINSIZE=$3
  CUTOFF=$4
  for i in $(ls $DATA/results_SIFT2/); do
    qsub -q standby -v SUBJ=${i},SESS=${SESS},WINSIZE=${WINSIZE},CUTOFF=${CUTOFF} get_nets.sub
  done
elif [ "$1" == "mod" ]; then
  SESS=$2
  WINSIZE=$3
  for i in $(ls $DATA/results_SIFT2/); do
    qsub -q standby -v SUBJ=${i},SESS=${SESS},WINSIZE=${WINSIZE} get_mod_trans.sub
  done
fi
