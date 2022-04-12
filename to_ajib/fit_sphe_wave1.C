#include "MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 12;
    Cal.fixZero = false;
    Cal.n_peaks = 8; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 = 0.01; // amplitude baseline
    Cal.mean0 = 1000; // average baseline
    Cal.sigma0 = 500; // sigma baseline
     
    Cal.peak1 = 0.025; // same for 1 p.e.
    Cal.mean1 = 6000;
    Cal.sigma1 = 1000;
    
    Cal.startpoint = 0.01; // amplitude for 2 p.e. 
        
    Cal.xmin = -1000; // range for graph display (not fit)
    Cal.xmax = 100000;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
