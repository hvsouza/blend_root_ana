#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 40;
    Cal.fixZero = false;
    Cal.n_peaks = 7; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.014; // amplitude baseline
    Cal.mean0 = -0.006; // average baseline
    Cal.sigma0 = 0.03; // sigma baseline
     
    Cal.peak1 = 0.014; // same for 1 p.e.
    Cal.mean1 = 0.13;
    Cal.sigma1 = 0.05;
    
    Cal.startpoint = 0.01; // amplitude for 2 p.e. 
        
    Cal.xmin = -0.7; // range for graph display (not fit)
    Cal.xmax = 5;

    Cal.deltaminus = 2.1;
    Cal.deltaplus = 1.3;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
