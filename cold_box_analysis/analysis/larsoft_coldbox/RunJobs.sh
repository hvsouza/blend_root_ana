#!/bin/bash

FCL_FILES=${FHICL_CB_CONFIG}
DATAPATH=/pnfs/dune/tape_backed/dunepro/vd-coldbox-top/full-reconstructed/2022/detector/test/VD_coldbox_TDE_2021/00/00/04/29
OUTPATH=/dune/data/users/ykermaid/coldbox/data/anatree
USERPATH=/dune/app/users/ykermaid/coldbox/data/anatree

cd ${DATAPATH}
FILES=( $(ls ./*.root | sort -n) ) 
cd ${USERPATH}

SEQUENCE=( $(seq 0 10 500) );

#source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh 
#setup jobsub_client

for (( i=0; i<${#SEQUENCE[@]}; i++ )); do 
	SUBFILES=( ${FILES[@]:${SEQUENCE[$i]}:10} )

        LIST=""
        for FILE in ${SUBFILES[@]}; do
            LIST+="${FILE}+"
        done

#        jobsub_submit -G dune -M -N 1 --memory=1800MB --disk=2GB --expected-lifetime=1h --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=${USERPATH}/RunAnaTree_$i.tar.gz -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105)' file://${USERPATH}/RunAnaTree.sh ${DATAPATH} ${OUTPATH} ${USERPATH} ${LIST} ${FCL_FILES}
        bash ${USERPATH}/RunAnaTree.sh ${DATAPATH} ${OUTPATH} ${USERPATH} ${LIST} ${FCL_FILES} 

done
