#!/usr/bin/env sh
# ________________________________________ #
# Author: Henrique Souza
# Filename: samples.C
# Created: 2021
# ________________________________________ #
#
# insert this line to .bashrc
# alias myclass="source ~__USER_PATH__/load_my_class.sh"
# to use the sample (analyzer) class, just type:
# myclass
# If you have different numbers of pts per waveform or file name, you can use, for instance:
# myclass 1000 "folder/file.root"

if [ -z $2 ]; then
    npts=5000
else
    npts=$2
fi

if [ -z $1 ]; then
    file='analyzed.root'
else
    file=$1
fi

if [ $file == 'no' ]; then
    eval 'root -e "#define memorydepth $npts" -e ".L __USER_PATH__/MYCODES.h"'
else
    eval 'root -e "#define memorydepth $npts" -e ".L __USER_PATH__/MYCODES.h" -e "ANALYZER s(\"s\")" -e "s.setAnalyzer(\"$file\")"'
fi