#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void fourier_tests(){
  WIENER wf("wn");
  WIENER wled("led");
  WIENER wpe("pe");
  WIENER wres("res");

  
  //A function to sample
  TF1 *fsin = new TF1("fsin", "sin(2*TMath::Pi()*10*x/(1e3))", 0, memorydepth*wf.step);
  fsin->SetNpx(10000);
  TCanvas *cprov = new TCanvas();
  fsin->Draw();
  
  Double_t x;
  
  wf.frequency = 250;
  wf.step = 4;
  
  TH1D* hnoise = new TH1D("hnoise","hnoise",memorydepth,0,memorydepth*wf.step);
  for(Int_t i = 0; i< memorydepth; i++){
      // hnoise->SetBinContent(i+1, fsin->Eval(i*4));   
      hnoise->SetBinContent(i+1,gRandom->Gaus(0,20));
  }
  wf.fft(hnoise);
  // cout << wf.hfft->Integral() << " " << wf.hfft->Integral("width") << " " << wf.hfft->GetBinWidth(1) << " " << endl;


  TCanvas *cwn = new TCanvas();
  wf.hfft->Draw();
  
  // TFile *f = new TFile("led_response.root","READ");
  // TH1D *hled = (TH1D*)f->Get("averaged_led");

  // TFile *f2 = new TFile("few_pe_response.root","READ");
  // TH1D *hpe = (TH1D*)f2->Get("averaged_few_pe");

  
  // wled.fft(hled);
  // wpe.fft(hpe);

  // TCanvas *cfft = new TCanvas();
  // wled.hfft->Draw();
  // wpe.hfft->Draw("SAME");


  // wres.deconvolve(wled,wpe,90);
  // // wres.shift_waveform(wres.hwvf,1200/4);

  // TCanvas *cres = new TCanvas();
  // wres.hwvf->Draw("");

  // cfft->cd();
  // wres.hfft->Draw("SAME");


  



  
}
