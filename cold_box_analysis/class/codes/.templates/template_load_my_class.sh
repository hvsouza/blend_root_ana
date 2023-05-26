#!/usr/bin/env sh
# ________________________________________ #
# Author: Henrique Souza
# Filename: samples.C
# Created: 2021
# ________________________________________ #
#
# insert this line to .bashrc
# alias myclass="source __USER_PATH__/load_my_class.sh"
# to use the sample (analyzer) class, go to any folder which have "analyzed.root" file and run:
# myclass
# If you want to execute a specific file instead, you can run
# myclass {path_to_file}/your_file.root 
# If you want to change old data, you need to give the number of points per waveform as following:
# myclass {file.root} {# of pts}
# If you want just to load the class without any data, you can run
# myclass no

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
