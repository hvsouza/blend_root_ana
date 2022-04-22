#!/bin/bash

FILEPATH=$1
OUTPATH=$2
USERPATH=$3
FILE=$4

NEWFORMAT=${FILE/.hdf5}.root
NEWFORMAT2=${FILE/.hdf5}_reco.root
NEWFORMAT3=${FILE/.hdf5}_reco_tree.root

cd $USERPATH

cp ${FILEPATH}/${FILE} ${OUTPATH}/
echo "lar -c vdcoldbox_raw_dataprep.fcl ${OUTPATH}/${FILE} -o ${OUTPATH}/${NEWFORMAT}"

lar -c vdcoldbox_raw_dataprep.fcl ${OUTPATH}/${FILE} -o ${OUTPATH}/${NEWFORMAT}

echo "lar -c standard_reco_vdcb_48deg_bot.fcl ${OUTPATH}/${NEWFORMAT} -o ${OUTPATH}/${NEWFORMAT2}"

lar -c standard_reco_vdcb_48deg_bot.fcl ${OUTPATH}/${NEWFORMAT} -o ${OUTPATH}/${NEWFORMAT2}

echo "lar -c standard_anatree_vdcb1_data.fcl ${OUTPATH}/${NEWFORMAT2} -o ${OUTPATH}/${NEWFORMAT3}"

lar -c standard_anatree_vdcb1_data.fcl ${OUTPATH}/${NEWFORMAT2} -o ${OUTPATH}/${NEWFORMAT3}

echo "mv ${USERPATH}/ana_hist.root ${OUTPATH}/${NEWFORMAT3}"
mv ${USERPATH}/ana_hist.root ${OUTPATH}/${NEWFORMAT3}
echo ""
