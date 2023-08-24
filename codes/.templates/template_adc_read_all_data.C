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
    
    r.baselineTime = 20000; // time limit for baseline
    r.chargeTime = 10300; // last time to integrate
    r.startCharge = 9950;
    r.maxRange = 5500; // max range to search for amplitude peak
    r.fast = 200; // fprompt fast integration time
    r.slow = 1700; //fprompt slow integration time
    r.exclusion_baseline = 35; // filtered waveform, anything above here will do +exclusion window while evaluating the baseline
    // r.exclusion_baselines = {35,15};
    r.exclusion_window = 1000; // time in ns that it will jump for baseline
    r.filter = 14; // denoise filter.. if filter = 0, no denoise is done.
    r.OnlyOneEvent = false; // Do you want only one event? Choose it wisely (you can set stopEvent)
    r.stopEvent = 1000; //
    r.noBaseline = false; // if you dont want to calculate baseline
    // r.saveFilter = true;

    // r.channels = {1};
    // r.channels = {1,2,3,4};
      

    
    r.adc_read_all_data();

    
}

