#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void showFilters(){

  ANALYZER z;
  ANALYZER zs("sub");
  ANALYZER zb("bandpass");
  ANALYZER zp("bandcut");
  ANALYZER zf("sub_gaus");

  Double_t lowFreqcut = 0.4;
  Double_t highFreqcut = 5;
  Double_t denoise = 24;

  z.setAnalyzer();
  zs.setAnalyzer();
  zb.setAnalyzer();
  zp.setAnalyzer();
  zf.setAnalyzer();

  Int_t &kch = z.kch;

  // for (Int_t i = 0; i < z.t1->GetEntries(); i++) {
  for (Int_t i = 0; i < 1; i++) {
    z.getWaveform(i);
    zs.getWaveform(i);
    zb.getWaveform(i);
    zp.getWaveform(i);
    zf.getWaveform(i);
  }

  z.drawGraph("AL");
  z.gwvf->SetLineColor(kGreen);

  // ____________________ Subtraction filter
  zs.applyFreqFilter(lowFreqcut,"low");
  zs.drawGraph("SAME");
  zs.gwvf->SetLineColor(kRed);
  zs.gwvf->SetLineWidth(3);
  for (Int_t j = 0; j < memorydepth; j++) {
    z.ch[kch].wvf[j] = z.ch[kch].wvf[j] - zs.ch[kch].wvf[j];
  }
  z.applyDenoise(denoise);
  z.drawGraph("SAME");
  z.gwvf->SetLineWidth(3);


  // ____________________ Subtraction filter + Gaus
  z.getWaveform(0);
  TCanvas *c2 = new TCanvas("c2");
  z.drawGraph();
  z.gwvf->SetLineColor(kGreen);
  zf.applyFreqFilter(lowFreqcut,"gaus");
  zf.drawGraph("SAME");
  zf.gwvf->SetLineColor(kRed);
  zf.gwvf->SetLineWidth(3);
  for (Int_t j = 0; j < memorydepth; j++) {
    z.ch[kch].wvf[j] = z.ch[kch].wvf[j] - zf.ch[kch].wvf[j];
  }
  z.applyDenoise(denoise);
  z.drawGraph("SAME");
  z.gwvf->SetLineWidth(3);


  // ____________________ Bandpass
  z.getWaveform(0);
  TCanvas *c3 = new TCanvas("c3");
  z.drawGraph();
  z.gwvf->SetLineColor(kGreen);
  zb.applyFreqFilter(1,"high");
  zb.applyFreqFilter(highFreqcut,"low");
  zb.applyDenoise(16);
  zb.drawGraph("SAME");
  zb.gwvf->SetLineWidth(3);
  c3->BuildLegend();
  // zp.w->hfft->Draw();

  // ____________________ Bandcut
  z.getWaveform(0);
  TCanvas *c4 = new TCanvas("c4");
  z.drawGraph();
  z.gwvf->SetLineColor(kGreen);
  zp.applyBandCut(0.2, 0.8);
  zp.applyDenoise(denoise);
  zp.drawGraph("SAME");
  zp.gwvf->SetLineWidth(3);
  c4->BuildLegend();
  // zp.w->hfft->Draw();
}
