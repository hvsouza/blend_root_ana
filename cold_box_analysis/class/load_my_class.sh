#!/usr/bin/env sh
# ________________________________________ #
# Author: Henrique Souza
# Filename: samples.C
# Created: 2021
# ________________________________________ #
#
# insert this line to .bashrc
# alias myclass="source ~/Dropbox/APC_Paris/Root/cold_box_analysis/class/load_my_class.sh"
# to use the sample (analyzer) class, just type:
# myclass
# If you have different numbers of pts per waveform or file name, you can use, for instance:
# myclass 1000 "folder/file.root"

if [ -z $1 ]; then
    npts=5000
else
    npts=$1
fi

if [ -z $2 ]; then
    file='analyzed.root'
else
    file=$2
fi

if [ $file == 'no' ]; then
    eval 'root -e "#define memorydepth $npts" -e ".L /home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"'
else
    eval 'root -e "#define memorydepth $npts" -e ".L /home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h" -e "SAMPLE s(\"s\",\"$file\")"'
fi
