#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void adc_read_all_data_1252(){
    
    // system("rm files.txt");
    // system("ls -1 *.txt| sed -e 's/.txt$//' > files.txt");
    
    Read r;

    r.dtime = 4;
    r.nbits = 14;
    r.isBinary = false;
    
    r.baselineTime = 20000; // time limit for baseline
    r.chargeTime = 1700; // last time to integrate 
    r.startCharge = 900;
    r.maxRange = 1700; // max range to search for amplitude peak
    r.fast = 200; // fprompt fast integration time
    r.slow = 1700; //fprompt slow integration time
    r.exclusion_baseline = 0.005; // filtered waveform, anything above here will do +exclusion window
    r.exclusion_window = 500; // time in ns that it will jump for baseline
    r.filter = 0.002; // denoise filter.. if filter = 0, no denoise is done. 
    r.OnlyOneEvent = false; // Do you want only one event? Choose it wisely 
    r.stopEvent = 2;
      
    r.channels = {1,2};
    
    
    r.adc_read_all_data();

    
}

