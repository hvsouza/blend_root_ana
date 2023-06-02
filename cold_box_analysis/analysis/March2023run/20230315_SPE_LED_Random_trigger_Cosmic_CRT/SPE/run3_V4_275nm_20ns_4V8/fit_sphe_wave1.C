#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_sphe_wave1(Int_t ch = 4){
    
    Calibration Cal;
    
    string histogram = "analyzed_" + to_string(ch);
    Cal.rootFile = "sphe_histograms_Ch"+to_string(ch)+".root";

    Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)

    Cal.rebin = 100;
    // Cal.fixZero = false;
    Cal.make_free_stddevs = true; // starts with false, if good fitting, change to true
    Cal.searchParameters(histogram.c_str(), 2, true); // give a first search in the parameters.
    // Cal.n_peaks = 12; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
    // Cal.peak0 =  0.014; // amplitude baseline
    // Cal.mean0 = 0; // average baseline
    // Cal.sigma0 = 1300; // sigma baseline

    // Cal.peak1 = 0.016; // same for 1 p.e.
    // Cal.mean1 = 6000;
    // Cal.sigma1 = 1700;

    // Cal.startpoint = 0.015; // amplitude for 2 p.e.

    // Cal.xmin = -10000; // range for graph display (not fit)
    // Cal.xmax = 60000;

    Cal.deltaminus = 0;
    Cal.deltaplus = 1;

    Cal.fit_sphe_wave(histogram.c_str(),false);
    Cal.hcharge->SetName("V4 (DCem1.2) Ch1");
}
