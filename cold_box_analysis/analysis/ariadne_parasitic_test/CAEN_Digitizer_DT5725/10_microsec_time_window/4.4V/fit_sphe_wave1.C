#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 100;
    Cal.fixZero = false;
    Cal.n_peaks = 8; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.028; // amplitude baseline
    Cal.mean0 = -500; // average baseline
    Cal.sigma0 = 350; // sigma baseline
     
    Cal.peak1 = 0.022; // same for 1 p.e.
    Cal.mean1 = 1600;
    Cal.sigma1 = 500;
    
    Cal.startpoint = 0.013; // amplitude for 2 p.e. 
        
    Cal.xmin = -1000; // range for graph display (not fit)
    Cal.xmax = 20000;

    Cal.deltaplus = 1.1;
    Cal.deltaminus = 2;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch3.root";
    Cal.fit_sphe_wave("analyzed_3");  
}
