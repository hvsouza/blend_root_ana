#!/usr/bin/env sh
if [ -z $1 ]; then
    npts=5000
else
    npts=$1
fi
eval 'root -e "#define memorydepth $npts" -e ".L /home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h" -e "SAMPLE s"'
