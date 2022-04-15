#!/bin/bash

DATAPATH=$1
OUTPATH=$2
USERPATH=$3
SUBFILES=$4
FCL_FILES=$5

SUBFILES=${SUBFILES//+/ }

source ${USERPATH}/.larsoft_profile

for FILE in ${SUBFILES[@]}; do

  echo "lar -c ${FCL_FILES}/standard_anatree_vdcb1_data.fcl ${DATAPATH}/${FILE} -o ${OUTPATH}/${FILE}"
  lar -c ${FCL_FILES}/standard_anatree_vdcb1_data.fcl ${DATAPATH}/${FILE} -o ${OUTPATH}/${FILE}

  echo "mv ${USERPATH}/ana_hist.root ${OUTPATH}/${FILE}"
  mv ${USERPATH}/ana_hist.root ${OUTPATH}/${FILE}
  echo ""

done
