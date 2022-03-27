#!/bin/bash

rsync -avzuh$1 --progress --exclude=*.pdf --exclude=*.zip --exclude=*.dat --exclude=*.txt --exclude=*.root --exclude=analyzed.root --exclude=sphe_waveforms_Ch*.root --exclude=sphe_histograms_darkCount_Ch* --exclude=*.so --exclude=*.d --exclude=*.pcm --exclude=*.trc --exclude=*.tar.gz *  "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/analysis"


rsync -avzuh$1 --progress --include="*/" --include="*graphs*/**" --exclude="*" --prune-empty-dirs . "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/analysis"
