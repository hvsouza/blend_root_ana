// ________________________________________ //
// Author: Henrique Souza
// Filename: giveMeSphe_darkCount.C
// Created: 2021
// ________________________________________ //
#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"



void giveMeSphe_darkCount(){
    
    SPHE dark;
    
    dark.led_calibration = true;
    
    dark.just_a_test = false;
    dark.just_this = 200;
    
    dark.tolerance = 5; // n sigmas (smoothed)
    dark.baseLimit = 3; // higher then this wont contribute to the baseline abs(baseLimit)
    dark.baselineTime = 5000;

    dark.start = 9950; // start the search for peaks or start the integration (led)
    dark.finish = 10250; // fisish the search or finish the integration (led)

    dark.timeLimit = 0; // time after LED signal
    dark.timeLow = 60; // integration time before peak
    dark.timeHigh = 400; // integration time after peak
    
    dark.lowerThreshold = -1; // threshold to detect noise (normal waveform)
    dark.maxHits = 1; // maximum hit before discarding  
    
    dark.too_big = 1000; // if there is a peak > "too_big" .. wait "waiting" ns
    dark.waiting = 1000;
    
    dark.filter = 14;
    dark.interactions = 45; // for moving avarage
    dark.shifter = 20;

    dark.dtime = 4.;

    dark.get_wave_form = false;
    dark.mean_before = 120; // time recorded before and after the peak found 
    dark.mean_after = 1000;
    dark.sphe_charge_ch0 = 1809.52; // wave0
    dark.sphe_charge2_ch0 = 3425.95; // wave0

    dark.channel = 1;
        
    dark.giveMeSphe_darkCount("analyzed");

    
    gROOT->SetBatch(kFALSE);
    
}
