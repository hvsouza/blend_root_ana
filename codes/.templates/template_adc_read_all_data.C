// ________________________________________ //
// Author: Henrique Souza
// Filename: adc_read_all_data.C
// Created: 2021
// ________________________________________ //
#define memorydepth 5000
#include "__USER_PATH__/MYCODES.h"


void adc_read_all_data(string datadir = "./"){
    
    system("rm files.log");
    int datadirlength = datadir.length();
    if(datadir[datadirlength-1] != '/'){
        datadir = datadir + "/";
    }
    string packdata = Form("ls -1 %s*.dat | sed -e 's/.dat$//' > files.log", datadir.c_str());
    // string packadata = Form("ls -1 %s*.dat | grep wave | sed -e 's/.dat$//' > files.log", datadir.c_str()); // to only check one channel
    // string packadata = Form("ls -1 %s*.dat | grep wave[2,3] | sed -e 's/.dat$//' > files.log", datadir.c_str()); // to only check channels 2 and 3
    system(packdata.c_str());
    
    Read r;
        
    r.dtime = 4; // step time in ns
    r.nbits = 14; // this is used only for baseline histogram, changing it to 12 might help
    r.isBinary = true;
    
    r.baselineTime = 1200; // time limit for baseline
    r.chargeTime = 3800; // last time to integrate
    r.startCharge = 1300;
    r.maxRange = 1600; // max range to search for amplitude peak
    r.fast = 100; // fprompt fast integration time
    r.slow = 5000; //fprompt slow integration time
    r.exclusion_baselines = {20, 20}; // filtered waveform, anything above here will trigger +exclusion window (two channels are being set to the same value)
    r.exclusion_window = 100; // time in ns that it will jump for baseline
    r.filter = 4; // denoise filter.. if filter = 0, no denoise is done. r.OnlyOneEvent = false; // Do you want only one event? Choose it wisely
    // r.OnlyOneEvent = true;
    r.stopEvent = 1000; // in case you want to collect only n events, set `stopEvent` = n and set `OnlyOneEvent` = True
    r.noBaseline = false; // Skip baseline correction when set to true
    r.polarity = -1; // Pulse polarity 1: positive, -1: negative
    // r.saveFilter = true; // This will save the waveforms with filter (not recommended)

    
    r.adc_read_all_data();

    
}

