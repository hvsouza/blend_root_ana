#!/bin/bash

# Quick executable for WaveDump, created by Henrique Souza
# This executable will create a folder the current date as name and execute the wavedump program inside it

# WARNING: the WaveDump program will always overwrite data files. Try to change your file name before running one acquisition

cd ~/home/lableptons/Desktop/WaveDumpData/

#now=$(date '+%Y%m%d %k:%M:%S')

#mkdir "$now"

#cd "$now"

gedit /etc/wavedump/WaveDumpConfig.txt &

wavedump
