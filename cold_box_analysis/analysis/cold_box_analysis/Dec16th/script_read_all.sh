#!/bin/bash

for x in $(eval "ls -d */ | grep Arapuca | sed 's:/$::'"); do
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
  if [ $variables -eq 938 ]; then
    eval "root -l -b -q ../adc_read_all_data_938.C"
  else
    eval "root -l -b -q ../adc_read_all_data_11876.C"
  fi
  cd "../"
done 
