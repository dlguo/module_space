#!/bin/sh

SESS=$1
WINSIZE=$2
CUTOFF=$3
for i in $(ls $DATA/results_SIFT2/); do
  qsub -q standby -v SUBJ=${i},SESS=${SESS},WINSIZE=${WINSIZE},CUTOFF=${CUTOFF} get_nets.sub
done

# SESS=$1
# WINSIZE=$2
# for i in $(ls $DATA/results_SIFT2/); do
#   qsub -q standby -v SUBJ=${i},SESS=${SESS},WINSIZE=${WINSIZE} get_mod_trans.sub
# done
