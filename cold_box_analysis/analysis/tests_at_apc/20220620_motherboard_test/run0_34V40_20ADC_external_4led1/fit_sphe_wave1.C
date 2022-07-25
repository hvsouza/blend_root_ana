#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 200;
    Cal.fixZero = false;
    Cal.n_peaks = 16; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.005; // amplitude baseline
    Cal.mean0 = -1; // average baseline
    Cal.sigma0 = 30; // sigma baseline
     
    Cal.peak1 = 0.008; // same for 1 p.e.
    Cal.mean1 = 200;
    Cal.sigma1 = 20;
    
    Cal.startpoint = 0.012; // amplitude for 2 p.e. 
        
    Cal.xmin = -800; // range for graph display (not fit)
    Cal.xmax = 6000;

    Cal.deltaplus = 1.2;
    Cal.deltaminus = 2;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
