#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fit_italy(){
    
  Calibration Cal;

  Cal.rebin = 6;
  Cal.make_free_stddevs = true;
  Cal.rootFile = "sphe_histograms_italy.root";
  Cal.searchParameters("wave_35V70_200ADC_1",2,true);

  // Cal.fixZero = false;
  // Cal.n_peaks = 8; // n peaks after baseline and 1 p.e.. Total = 2 + n_peaks
  // Cal.peak0 =  0.014; // amplitude baseline
  // Cal.mean0 = 0; // average baseline
  // Cal.sigma0 = 1300; // sigma baseline

  // Cal.peak1 = 0.016; // same for 1 p.e.
  // Cal.mean1 = 6000;
  // Cal.sigma1 = 1700;

  // Cal.startpoint = 0.015; // amplitude for 2 p.e.

  // Cal.xmin = -10000; // range for graph display (not fit)
  // Cal.xmax = 60000;

  // Cal.deltaplus = 1.35;
  // Cal.deltaminus = 1.4;
  Cal.fit_sphe_wave("wave_35V70_200ADC_1",false);
}
