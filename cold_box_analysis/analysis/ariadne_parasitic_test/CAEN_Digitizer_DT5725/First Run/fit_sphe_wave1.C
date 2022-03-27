#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"
void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 80;
    Cal.fixZero = true;
    Cal.n_peaks = 4; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 = 0; // amplitude baseline
    Cal.mean0 = 0; // average baseline
    Cal.sigma0 = 0; // sigma baseline
     
    Cal.peak1 = 0.02; // same for 1 p.e.
    Cal.mean1 = 2100;
    Cal.sigma1 = 400;
    
    Cal.startpoint = 0.018; // amplitude for 2 p.e. 
        
    Cal.xmin = 1200; // range for graph display (not fit)
    Cal.xmax = 100000;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch3.root";
    Cal.fit_sphe_wave("analyzed_3");  
}
