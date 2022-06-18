#!/bin/bash

read -p "Are you sure (y/n)? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    return
fi

anadir=$(pwd)
datadir=~/Documents/ADC_data/coldbox_data
cd $datadir
for x in $(eval "ls -d  */ | grep 20220614"); do
    cd $x;
    for subruns in $(eval "ls -d  */ | grep run"); do
        # echo $subruns
        cd $subruns
        eval "root -l -b -q $anadir/adc_read_all_data.C"
        mkdir -p $anadir/$x$subruns
        mv analyzed.root $anadir/$x$subruns
        cd ../
    done
done

cd $anadir
