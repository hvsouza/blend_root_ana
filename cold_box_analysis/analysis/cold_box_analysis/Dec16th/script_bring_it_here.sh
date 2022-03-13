#!/bin/bash
 
for x in $(eval "ls -1 *.txt| sed -e 's/.txt$//'"); do
  eval "mkdir $x"
  eval "mv $x.txt $x"
done 
  
