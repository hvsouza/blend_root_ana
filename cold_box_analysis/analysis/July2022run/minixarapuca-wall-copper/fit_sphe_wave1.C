#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 2; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 50;
    Cal.fixZero = false;
    Cal.n_peaks = 12; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.004; // amplitude baseline
    Cal.mean0 = -200; // average baseline
    Cal.sigma0 = 1000; // sigma baseline
     
    Cal.peak1 = 0.007; // same for 1 p.e.
    Cal.mean1 = 2000;
    Cal.sigma1 = 1000;
    
    Cal.startpoint = 0.009; // amplitude for 2 p.e. 
        
    Cal.xmin = -1500; // range for graph display (not fit)
    Cal.xmax = 20000;

    Cal.deltaplus = 1.5;
    Cal.deltaminus = 1.5;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
