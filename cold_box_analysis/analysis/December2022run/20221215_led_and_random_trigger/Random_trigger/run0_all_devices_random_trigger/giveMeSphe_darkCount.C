#define memorydepth 250000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"



void giveMeSphe_darkCount(){
    
    SPHE dark;
    
    dark.led_calibration=false;
    
    dark.just_a_test = false;
    dark.just_this = 200;
    // dark.check_selection = false;

    dark.fix_threshold = true;
    dark.tolerance = 7.3; // n sigmas (smoothed)
    dark.baseLimit = 18; // higher then this wont contribute to the baseline abs(baseLimit)
    dark.baselineTime = 500000;
    dark.social_distance = 1;
    // dark.start = 10310; // start the search for peaks or start the integration (led)
    // dark.finish = 10700; // fisish the search or finish the integration (led)

    dark.start = 0; // start the search for peaks or start the integration (led)
    dark.finish = 1000000; // fisish the search or finish the integration (led)

    dark.timeLimit = 0; // time after LED signal
    dark.timeLow = 8; // integration time before peak
    dark.timeHigh = 412; // integration time after peak
    
    dark.lowerThreshold = -7; // threshold to detect noise (normal waveform)
    dark.maxHits = 5; // maximum hit before discarding
    
    dark.too_big = 300; // if there is a peak > "too_big" .. wait "waiting" ns
    dark.waiting = 1000;

    dark.filter = 16;
    dark.interactions = 45; // for moving avarage
    dark.shifter = 20;
    
    dark.dtime = 4.;

    dark.get_wave_form = false;
    dark.mean_before = 120; // time recorded before and after the peak found 
    dark.mean_after = 1000;
    // dark.sphe_charge_ch0 = 4486.84; //wave0
    // dark.sphe_charge2_ch0 = 8910.46; // wave0
    dark.deltaplus = 1.3;
    dark.deltaminus = 1.5;
    

    dark.sphe_charge_ch0 = 5343.7; //wave0
    dark.sphe_charge2_ch0 = 10665.5; // wave0

    dark.channel = 5;
    dark.check_selection = false;
        
    
    dark.giveMeSphe_darkCount("analyzed");


    
    gROOT->SetBatch(kFALSE);
    
}
