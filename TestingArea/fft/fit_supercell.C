#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"
// #include "/home/henrique/Dropbox/Unicamp/Doutorado/Root/Programs/italy/ADC_LAr_SuperCell/class/MYCODES.h"
void fit_supercell(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 50;
    Cal.make_free_stddevs = true;
    // Cal.fixZero = false;
    // Cal.n_peaks = 8; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    // Cal.peak0 = 0.01; // amplitude baseline
    // Cal.mean0 = 1000; // average baseline
    // Cal.sigma0 = 500; // sigma baseline

    // Cal.peak1 = 0.025; // same for 1 p.e.
    // Cal.mean1 = 6000;
    // Cal.sigma1 = 1000;

    // Cal.startpoint = 0.01; // amplitude for 2 p.e.

    // Cal.xmin = -1000; // range for graph display (not fit)
    // Cal.xmax = 100000;

    Cal.rootFile = "sphe_histograms_darkCount_Ch1.root";
    Cal.searchParameters("analyzed_1",2,true);
    Cal.fit_sphe_wave("analyzed_1");  
}
