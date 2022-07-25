#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void persistance_plots(){
  const Int_t n = 6;
  Int_t channels[n] = {1,2,3,4,5,6};
  // TFile *fin = new TFile("cathodexarapuca/analyzed.root","READ");
  // TFile *fin = new TFile("minixarapuca-wall-copper/analyzed.root","READ");
  TFile *fin = new TFile("pulser/analyzed.root","READ");
  TTree *t1 = (TTree*)fin->Get("t1");

  vector<string> names = {"miniArapuca (A1ch1)","miniArapuca (A1ch2)","miniArapuca (Good ch1)","xA2-ch2 (DCem46V)","xA2-ch1 (DCem46V)","xA1-ch2 (DCem36V)"};
  vector<TCanvas *> c(n);

  vector<TH2D*> h(n);

  for(Int_t i = 0; i<6; i++){
    c[i] = new TCanvas(Form("c%i",i+1),Form("c%i",i+1));
    c[i]->cd();
    h[i] = new TH2D(Form("h%i",i+1),names[i].c_str(),5000,0,10000,5000,-20,700);
    t1->Draw(Form("Ch%i.wvf[]*2000/pow(2,14):Iteration$*2>>h%i",channels[i],channels[i]),"","colz");
    h[i]->GetXaxis()->SetRangeUser(3500,7000);
    h[i]->GetXaxis()->SetTitle("Time (ns)");
    h[i]->GetYaxis()->SetTitle("Amplitude (mV)");
    h[i]->SetLineWidth(4);
    h[i]->SetStats(kFALSE);
    c[i]->BuildLegend();
  }
  
}
