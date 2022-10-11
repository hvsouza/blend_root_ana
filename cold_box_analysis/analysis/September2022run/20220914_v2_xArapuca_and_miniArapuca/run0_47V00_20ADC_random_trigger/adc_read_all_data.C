#define memorydepth 250000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void adc_read_all_data(){
    
    system("rm files.log");
    system("ls -1 *.dat| sed -e 's/.dat$//' > files.log");
    
    Read r;
        
    r.dtime = 4;
    r.basebits = 14;
    r.isBinary = true;
    
    r.baselineTime = 250000*4; // time limit for baseline
    r.chargeTime = 5500; // last time to integrate
    r.startCharge = 5320;
    r.maxRange = 5500; // max range to search for amplitude peak
    r.fast = 200; // fprompt fast integration time
    r.slow = 1700; //fprompt slow integration time
    r.exclusion_baseline = 35; // filtered waveform, anything above here will do +exclusion window
    r.exclusion_window = 1000; // time in ns that it will jump for baseline
    r.filter = 14; // denoise filter.. if filter = 0, no denoise is done.
    r.OnlyOneEvent = false; // Do you want only one event? Choose it wisely
    r.stopEvent = 1000;
    r.noBaseline = false;
    // r.saveFilter = true;
      
    // r.channels = {5,6};
    r.channels = {1,2,5,6};
    
    r.adc_read_all_data();

    
}

