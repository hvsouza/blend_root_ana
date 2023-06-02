#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void adc_read_all_data_led(){
    
    system("rm files.log");
    // system("ls -1 *.dat| sed -e 's/.dat$//' > files.log");

    system("ls -1 *.dat | grep wave[3,4] | sed -e 's/.dat$//' > files.log"); // to only check one channel
    // system("ls -1 *.dat | grep wave[5,6] | sed -e 's/.dat$//' > files.log"); // to only check one channel
    Read r;
        
    r.dtime = 4;
    r.nbits = 14;
    r.isBinary = true;
    
    r.baselineTime = 20000; // time limit for baseline
    r.chargeTime = 10500; // last time to integrate
    r.startCharge = 10100;
    r.maxRange = 10500; // max range to search for amplitude peak
    r.fast = 200; // fprompt fast integration time
    r.slow = 1700; //fprompt slow integration time
    r.exclusion_baselines = {35, 3, 20, 13, 17, 20, 4, 4}; // filtered waveform, anything above here will do +exclusion window
    r.exclusion_window = 1000; // time in ns that it will jump for baseline
    r.filter = 16; // denoise filter.. if filter = 0, no denoise is done. r.OnlyOneEvent = false; // Do you want only one event? Choose it wisely
    // r.OnlyOneEvent = true;
    r.stopEvent = 100;
    r.noBaseline = false;
    // r.saveFilter = true;
      

    r.adc_read_all_data();

    
}

