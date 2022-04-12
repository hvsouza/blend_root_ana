#!/bin/bash

read -p "Are you sure (y/n)? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    return
fi

for x in $(eval "ls -d */ | grep V"); do
  cd $x;
  eval "root -l -b -q ../adc_read_all_data.C"
  # rm *.dat
  cd "../"
done 
