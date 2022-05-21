#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fourier_test(){
  WIENER wmu("mu");
  WIENER wpe("pe");
  WIENER wres("res");
  
  TFile *fwvf = new TFile("ch2_summed_waveform.root","READ");
  TFile *fsphe = new TFile("Merged_LED_ch2_3k_3p5k.root","READ");

  TH1D *hmu = (TH1D*)fwvf->Get("total_waveform");
  TH1D *hpe = (TH1D*)fsphe->Get("summed_wf");
  
  
  wmu.fft(hmu);
  wpe.fft(hpe);

  TCanvas *cfft = new TCanvas();
  wmu.hfft->Draw();
  wpe.hfft->Draw("SAME");


  wres.deconvolve(wmu,wpe,50);
  // wres.shift_waveform(0,1200)

  TCanvas *cres = new TCanvas();
  wres.hwvf->Draw("hist");

  cfft->cd();
  wres.hfft->Draw("SAME");


  



  
}
