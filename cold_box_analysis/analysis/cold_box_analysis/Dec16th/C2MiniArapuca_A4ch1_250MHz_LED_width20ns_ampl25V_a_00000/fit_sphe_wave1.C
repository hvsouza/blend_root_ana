#include "/home/henrique/Dropbox/Unicamp/Doutorado/Root/Programs/italy/ADC_LAr_SuperCell/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 500;
    Cal.fixZero = false;
    Cal.n_peaks = 4; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.025; // amplitude baseline
    Cal.mean0 = -0.1; // average baseline
    Cal.sigma0 = 0.05; // sigma baseline
     
    Cal.peak1 = 0.025; // same for 1 p.e.
    Cal.mean1 = 0.3;
    Cal.sigma1 = 0.05;
    
    Cal.startpoint = 0.015; // amplitude for 2 p.e. 
        
    Cal.xmin = -0.7; // range for graph display (not fit)
    Cal.xmax = 5;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
