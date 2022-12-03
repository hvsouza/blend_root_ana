// ________________________________________ //
// Author: Henrique Souza
// Filename: processBaseline.C
// Created: 2022-12-02
// ________________________________________ //

#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

Double_t baseline(Double_t v[],Int_t &selection);

TH1D *hbase = new TH1D("hbase", "hbase",TMath::Power(2,14),0,TMath::Power(2,14));

Double_t dtime = 4;
Double_t baselineTime = memorydepth*dtime;
void processBaseline(){

  ANALYZER z("z");
  baselineTime = 10000;

  z.setAnalyzer();
  z.kch = 1;

  TH1D *h = new TH1D("h", "h",TMath::Power(2,14),0,TMath::Power(2,14));
  TH1D *hres = new TH1D("hres","hres",200,0,0);
  TH1D *hstd = new TH1D("hstd","hstd",200,0,50);
  TH1D *hcharge = new TH1D("hcharge", "hcharge",50000,0,0);

  for (Int_t i = 0; i < 1; i++) {
    z.getWaveform(i,z.kch);
    z.applyDenoise(16);
    h->Reset();
    for (Int_t j = 0; j < baselineTime/dtime; j++) {
      h->Fill(z.ch[z.kch].wvf[j]);
    }
    Double_t res0 = h->GetBinCenter(h->GetMaximumBin());
    hres->Fill(res0);
    hstd->Fill(h->GetStdDev());
    // for (Int_t j = 0; j < baselineTime/dtime; j++) {
    //   h->Fill(z.ch[z.kch].wvf[j]);
    // }

  }

  TCanvas *c1 = new TCanvas("c1","c1",1920,0,800,800);
  h->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",1920,0,800,800);
  hstd->Draw();
  // h->Draw();
  // Double_t res0 = h->GetBinCenter(h->GetMaximumBin());
  // TLine *l = new TLine(res0,0,res0,h->GetMaximum());
  // l->SetLineColor(kRed);
  // l->Draw("SAME");

}

Double_t baseline(Double_t v[],Int_t &selection, Int_t idx, Int_t mevent){
  // Double_t result = 0;
  // hbase->Reset();
  // for(Int_t i=0; i<baselineTime/dtime; i++) hbase->Fill(v[i]);
  // Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());
  // Double_t hmean = hbase->GetMean();
  // Double_t hstd = hbase->GetStdDev();


  // Bool_t changed_mean = false;

  // Double_t bins=0;
  // for(Int_t i=0; i<baselineTime/dtime;){
  //   if(v[i] > res0 + exclusion_baseline || v[i]<res0 - exclusion_baseline) {
  //     i+=exclusion_window/dtime;
  //   }
  //   else{
  //     result += v[i];
  //     bins+=1;
  //     i++;
  //   }
  // }
  // if(bins>0)result/=bins;
  // if(bins > (baselineTime/dtime)/3.){
  //   selection = 0;
  //   return result;
  // }
  // else{
  //   if(changed_mean==false)selection = 1;
  //   else selection = 2;
  //   return res0;
  // }
  return 0;
}
