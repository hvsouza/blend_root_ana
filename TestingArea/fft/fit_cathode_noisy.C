#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_cathode_noisy(){
    
    Calibration Cal;
    
    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    
    Cal.rebin = 125;
    Cal.make_free_stddevs = true;
    Cal.rootFile = "sphe_histograms_darkCount.root";
    Cal.searchParameters("analyzed_1",2,true);

    // Cal.fixZero = false;
    // Cal.n_peaks = 5; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    // Cal.peak0 =  0.024; // amplitude baseline
    // Cal.mean0 = 0; // average baseline
    // Cal.sigma0 = 400; // sigma baseline

    // Cal.peak1 = 0.006; // same for 1 p.e.
    // Cal.mean1 = 1000;
    // Cal.sigma1 = 400;

    // Cal.startpoint = 0.004; // amplitude for 2 p.e.

    Cal.xmin = -800; // range for graph display (not fit)
    // Cal.xmax = 60000;

    // Cal.deltaplus = 1.35;
    // Cal.deltaminus = 1.4;
    Cal.fit_sphe_wave("analyzed_1",false);
}
