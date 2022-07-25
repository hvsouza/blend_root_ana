#!/bin/bash

set -e # if something fails, crash everthing and dont keep going

read -p "Are you sure (y/n)? " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    return
fi

anadir=$(pwd)
datadir=~/Documents/ADC_data/coldbox_data/June2022run
cd $datadir
for x in $(eval "ls -d  */ | grep 20220614"); do
# for x in $(eval "ls -d  */ | grep 20220615_cosmic_muons_data_cathode_off"); do
# for x in $(eval "ls -d  */ | grep 20220615_cosmic_muons_data_cathode_on"); do
# for x in $(eval "ls -d  */ | grep 20220615_L"); do
    cd $x
    for subruns in $(eval "ls -d  */ | grep run"); do
        # echo $x$subruns
        cd $subruns
        echo $(pwd)
        eval "root -l -b -q $anadir/adc_read_all_data.C"
        mkdir -p $anadir/$x$subruns
        mv analyzed.root $anadir/$x$subruns
        cd ../
    done
    cd ../
done

cd $anadir

set +e # now back to normal
