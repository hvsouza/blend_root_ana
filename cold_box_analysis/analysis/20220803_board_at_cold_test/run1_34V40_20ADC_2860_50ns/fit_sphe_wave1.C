#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 125;
    Cal.fixZero = false;
    Cal.make_free_stddevs = true;

    Cal.n_peaks = 6; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    Cal.peak0 =  0.007; // amplitude baseline
    Cal.mean0 = 0; // average baseline
    Cal.sigma0 = 900; // sigma baseline
     
    Cal.peak1 = 0.015; // same for 1 p.e.
    Cal.mean1 = 3400;
    Cal.sigma1 = 1100;
    
    Cal.startpoint = 0.01; // amplitude for 2 p.e.
        
    Cal.xmin = -5000; // range for graph display (not fit)
    Cal.xmax = 30000;

    Cal.deltaplus = 1.35;
    Cal.deltaminus = 1.4;
    
    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.fit_sphe_wave("analyzed_1");  
}
