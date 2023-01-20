#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"



void giveMeSphe_darkCount(){
    
    SPHE dark;
    
    dark.led_calibration=true;
    
    dark.just_a_test = false;
    dark.just_this = 200;
    // dark.check_selection = false;
    
    dark.tolerance = 3.5; // n sigmas (smoothed)
    dark.baseLimit = 2; // higher then this wont contribute to the baseline abs(baseLimit)
    dark.baselineTime = 5000;
    dark.start = 10320; // start the search for peaks or start the integration (led)
    dark.finish = 10600; // fisish the search or finish the integration (led)

    // dark.start = 10320; // start the search for peaks or start the integration (led)
    // dark.finish = 10500; // fisish the search or finish the integration (led)

    dark.timeLimit = 0; // time after LED signal
    dark.timeLow = 8; // integration time before peak
    dark.timeHigh = 340; // integration time after peak
    
    dark.lowerThreshold = -9999; // threshold to detect noise (normal waveform)
    dark.maxHits = 1; // maximum hit before discarding  
    
    dark.too_big = 1000; // if there is a peak > "too_big" .. wait "waiting" ns
    dark.waiting = 1000;

    dark.filter = 16;
    dark.interactions = 45; // for moving avarage
    dark.shifter = 20;
    
    dark.dtime = 4.;

    dark.get_wave_form = true;
    dark.mean_before = 120; // time recorded before and after the peak found 
    dark.mean_after = 1000;
    // dark.sphe_charge_ch0 = 4486.84; //wave0
    // dark.sphe_charge2_ch0 = 8910.46; // wave0
    // dark.deltaplus = 1.35;
    // dark.deltaminus = 1.4;

    dark.deltaplus = 1.6;
    dark.deltaminus = 1.25;

    // dark.sphe_charge_ch0 = 6005.51; //wave0
    // dark.sphe_charge2_ch0 = 12007.6; // wave0
    dark.sphe_charge_ch0 = 2434.57; //wave0
    dark.sphe_charge2_ch0 = 4688.44; // wave0
    
    dark.channel = 6;
        
    
    dark.giveMeSphe_darkCount("analyzed");


    
    gROOT->SetBatch(kFALSE);
    
}
