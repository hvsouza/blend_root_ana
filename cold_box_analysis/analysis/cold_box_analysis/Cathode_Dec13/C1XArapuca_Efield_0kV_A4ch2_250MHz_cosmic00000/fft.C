#define memorydepth 1124
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fft(){
  WIENER wmu("mu");
  WIENER wpe("pe");
  WIENER wres("res");

 
 
  TFile *f = new TFile("mu_response.root","READ");
  TH1D *hmu = (TH1D*)f->Get("averaged_mu");

  TFile *f2 = new TFile("few_pe_response.root","READ");
  TH1D *hpe = (TH1D*)f2->Get("averaged_few_pe");

  
  wmu.fft(hmu);
  wpe.fft(hpe);

  TCanvas *cfft = new TCanvas();
  wmu.hfft->Draw();
  wpe.hfft->Draw("SAME");


  wres.deconvolve(wmu,wpe,50);
  // wres.shift_waveform(wres.hwvf,1200/4);

  TCanvas *cres = new TCanvas();
  wres.hwvf->Draw("");

  cfft->cd();
  wres.hfft->Draw("SAME");


  



  
}
