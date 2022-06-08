#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

#include "TH1D.h"
#include "TVirtualFFT.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"


 
void testing_fft()
{
   
  TFile *fsphe = new TFile("few_pe_response.root","READ");
  TFile *fsphe2 = new TFile("few_pe_response1000.root","READ");
  TH1D *hsphe = (TH1D*)fsphe->Get("averaged_few_pe");
  TH1D *hsphe2 = (TH1D*)fsphe2->Get("averaged_few_pe");

  TFile *fwvf = new TFile("led_response.root","READ");
  TH1D *hsignal = (TH1D*)fwvf->Get("averaged_led");

   
  WIENER wpe("pe",4,250,1e-9,1e6,2500);
  WIENER wpe2("pe2",4,250,1e-9,1e6,2500);
  WIENER wpe3("pe3",4,250,1e-9,1e6,2500);
  WIENER back("back",4,250,1e-9,1e6,2500);
  WIENER wsignal("signal",4,250,1e-9,1e6,2500);
  WIENER wdec("dec",4,250,1e-9,1e6,2500);


  wpe.fft(hsphe);
  wsignal.fft(hsignal);
  wpe2.fft(hsphe2);
  wpe3.fft(hsphe2);
  wpe2.rescale_histogram(2500);
  
  back.backfft(wpe3);


  TCanvas *cfft = new TCanvas();

  wpe.hfft->Draw("");
  wpe2.hfft->SetLineColor(kRed);
  wpe2.hfft->Draw("SAME");
  // wpe3.hfft->Draw("SAME");

  TCanvas *cback = new TCanvas();
  // back.hwvf->Draw();
  hsphe->SetLineColor(kRed);
  hsphe->Draw("");
  back.hwvf->Draw("same");
  TCanvas *gdec = new TCanvas();
  wdec.deconvolve(wsignal,wpe,30);
  wdec.hwvf->Draw();
}
 
