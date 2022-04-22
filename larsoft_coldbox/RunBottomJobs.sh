#!/bin/bash

DATAPATH=/pnfs/dune/tape_backed/dunepro/vd-coldbox-bottom/raw/2022/detector/test/None/00/01
OUTPATH=/dune/data/users/hsouza/larsoft_coldbox/work/BDE
USERPATH=/dune/app/users/hsouza/larsoft_coldbox/work/BDE

# I need 3 fors here, because there are subdirectoties
cd ${DATAPATH}
DIR1=( $(ls) )
for i in ${DIR1[@]}; do
    
    cd ${DATAPATH}/$i
    # pwd 
    DIR2=$(ls)
    for j in ${DIR2[@]}; do
	cd $j
        FILE=($(ls))
        FILEPATH=$(pwd)
        bash ${USERPATH}/RunAnaTreeBottom.sh ${FILEPATH} ${OUTPATH} ${USERPATH} ${FILE}
	# echo $FILE
	# echo $FILEPATH
    done
done

cd $USERPATH
