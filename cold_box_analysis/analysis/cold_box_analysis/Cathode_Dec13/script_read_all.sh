#!/bin/bash

for x in $(eval "ls -d */ | grep 000"); do
  cd $x;
  eval "root -l -b -q ../adc_read_all_data.C"
  cd "../"
done 
