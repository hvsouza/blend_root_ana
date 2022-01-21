#!/bin/bash

for x in $(eval "ls -d */ | grep 000 | sed 's:/$::'"); do
  cd $x;
  counter=0
  for variables in $(sed '2q;d' $x.txt); do
    let counter+=1
#     if [ $counter -eq 4 ]; then
#       echo $variables
#     fi
  done
  variables="${variables//[$'\t\r\n ']}"
#   echo $variables
  if [ $variables -eq 2502 ]; then
    eval "root -l -b -q ../adc_read_all_data_2502.C"
  else
    eval "root -l -b -q ../adc_read_all_data_1252.C"
  fi
  cd "../"
done 
