#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void simple_fft(){

  TFile *fsphe = new TFile("SPE_HPK_AW.root","READ");
  TH1D *hsphe = (TH1D*)fsphe->Get("Run3_ch0_ADC0_0V_ScintProfFirstSignalBin_17_");
   TFile *flaser = new TFile("LASER_HPK_AW.root","READ");
  TH1D *hlaser = (TH1D*)flaser->Get("Run11_ch0_ADC0_0V_ScintProfFirstSignalBin_11_");

  TFile *falpha = new TFile("ALPHA_HPK_AW.root","READ");
  TH1D *halpha = (TH1D*)falpha->Get("Run27_ch0_ADC0_0V_ScintProfFirstSignalBin_5_");

  WIENER wsphe("wsphe",4,250,1e-9,1e6,memorydepth);
  WIENER wlaser("wlaser",4,250,1e-9,1e6,memorydepth);
  WIENER walpha("walpha",4,250,1e-9,1e6,memorydepth);
  WIENER wdec("wdec",4,250,1e-9,1e6,memorydepth);

  wsphe.fft(hsphe);
  walpha.fft(halpha);
  wlaser.fft(hlaser);
  wdec.deconvolve(walpha,wsphe,2);
  wdec.hwvf->Draw();

  wsphe.shift_waveform(hsphe,3000/4);
  wlaser.shift_waveform(hlaser,3000/4);
  walpha.shift_waveform(halpha,3000/4);

  hsphe->Scale(1/hsphe->GetMaximum());
  hlaser->Scale(1/hlaser->GetMaximum());
  halpha->Scale(1/halpha->GetMaximum());

  hlaser->SetLineColor(kRed);
  hlaser->SetLineWidth(2);
  
  hsphe->SetLineColor(kBlue);
  hsphe->SetLineWidth(2);
  
  halpha->SetLineColor(kBlack);
  halpha->SetLineWidth(2);

  hlaser->SetTitle("LASER");
  halpha->SetTitle("ALPHA");
  hsphe->SetTitle("SPE");

  
  TCanvas *c2 = new TCanvas();
  
  hlaser->Draw("hist");
  halpha->Draw("hist SAME");
  hsphe->Draw("hist SAME");
}
