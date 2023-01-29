#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void adc_read_all_data_cosmic(){
    
    system("rm files.log");
    system("ls -1 *.dat | grep wave7 | sed -e 's/.dat$//' > files.log");
    
    Read r;
        
    r.dtime = 4.;
    r.basebits = 14;
    r.isBinary = true;
    
    r.baselineTime = 2600; // time limit for baseline
    r.chargeTime = 7000; // last time to integrate
    r.startCharge = 2800;
    r.maxRange = 5000; // max range to search for amplitude peak
    r.fast = 300; // fprompt fast integration time
    r.slow = 4200; //fprompt slow integration time
    // r.exclusion_baselines = {35,15}; // filtered waveform, anything above here will do +exclusion window
    r.exclusion_baselines = {30,7,35,15}; // filtered waveform, anything above here will do +exclusion window
    r.exclusion_window = 1000; // time in ns that it will jump for baseline
    r.filter = 4; // denoise filter.. if filter = 0, no denoise is done.
    // r.OnlyOneEvent = true; // Do you want only one event? Choose it wisely
    r.stopEvent = 1000;
    // r.noBaseline = true;
    // r.saveFilter = true;
      
    r.channels = {1};
    // r.channels = {1,2,5,6};
    
    r.adc_read_all_data();

    
}

