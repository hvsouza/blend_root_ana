#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void makeHistograms(){

  ANALYZER z;

  ANALYZER zs("sub");
  ANALYZER zbase("base");
  ANALYZER zb("bandpass");
  ANALYZER zp("bandcut");
  ANALYZER zf("sub_gaus");

  Double_t lowFreqcut = 0.5;
  Double_t highFreqcut = 5;
  Double_t denoise = 16;

  z.setAnalyzer();
  zbase.setAnalyzer();
  zs.setAnalyzer();
  zb.setAnalyzer();
  zp.setAnalyzer();
  zf.setAnalyzer();

  Int_t &kch = z.kch;
  vector<string> filters = {"Raw", "Subtraction (Low pass)","Subtraction (Gaus)", "Bandpass", "Bandcut"};
  Int_t nfilters = filters.size();
  vector<TH1D*> hc(nfilters);
  vector<TCanvas*> c(nfilters);
  vector<Double_t> charge(nfilters,0);

  WIENER *wtemp = new WIENER("wtemp",4,250,1e-9,1e6,5000);

  fstream f("output.dat", ios::out | ios::binary);
  for(Int_t i = 0; i < nfilters; i++){
    hc[i] = new TH1D(filters[i].c_str(),filters[i].c_str(),50000,0,0);
    c[i] = new TCanvas(Form("c%d",i));
  }

  zbase.setFreqFilter(lowFreqcut,"low");
  zf.setFreqFilter(lowFreqcut,"gaus");
  wtemp->setBandCut(0.2, 0.8, wtemp);
  TH2D *hp = new TH2D("hp","hp",5000,0,20000,1000,-200,800);
  for (Int_t i = 0; i < z.t1->GetEntries(); i++) {
  bool makebreak = false;
  // for (Int_t i = 0; i < 500; i++) {
    makebreak = false;
    // printf("\rEvaluating... %d",i);
    // fflush(stdout);
    z.getWaveform(i);
    zbase.getWaveform(i);
    zs.getWaveform(i);
    zb.getWaveform(i);
    zp.getWaveform(i);
    zf.getWaveform(i);

    // ____________________ Subtraction filter
    zbase.applyFreqFilter();
    for (Int_t j = 0; j < memorydepth; j++) {
      zs.ch[kch].wvf[j] = zs.ch[kch].wvf[j] - zbase.ch[kch].wvf[j];
      if(zs.ch[kch].wvf[j] > 200) makebreak = true;
    }
    if (makebreak) continue;
    // zs.applyDenoise(denoise);
    for (Int_t j = 0; j < memorydepth; j++) {
      f.write((char*)&zs.ch[kch].wvf[j],8);
      hp->Fill(j*4,zs.ch[kch].wvf[j]);
    }
    

    // // ____________________ Subtraction filter + Gaus
    // zbase.getWaveform(i);
    // zf.applyFreqFilter();
    // for (Int_t j = 0; j < memorydepth; j++) {
    //   zf.ch[kch].wvf[j] = zf.ch[kch].wvf[j] - zbase.ch[kch].wvf[j];
    // }
    // zf.applyDenoise(denoise);

    // // ____________________ Bandpass
    // zb.setFreqFilter(1,"high");
    // zb.applyFreqFilter();
    // zb.setFreqFilter(highFreqcut,"low");
    // zb.applyFreqFilter();
    // zb.applyDenoise(16);

    // ____________________ Bandcut
    // zp.applyBandCut(wtemp);
    // zp.applyDenoise(denoise);



    for(Int_t j = 10370/z.dtime; j < 10500/z.dtime; j++){
      charge[0] += z.ch[kch].wvf[j];
      charge[1] += zs.ch[kch].wvf[j];
      charge[2] += zf.ch[kch].wvf[j];
      charge[3] += zb.ch[kch].wvf[j];
      charge[4] += zp.ch[kch].wvf[j];
    }

    for(Int_t k = 0; k < nfilters; k++)
    {
      hc[k]->Fill(charge[k]*z.dtime);
      charge[k] = 0;
    }
  }
  for(Int_t k = 0; k < nfilters; k++)
  {
    c[k]->cd();
    hc[k]->Draw();
  }

  TCanvas *cp = new TCanvas("cp");
  hp->Draw("colz");


  TFile *fout = new TFile("sphe_histograms_darkCount.root","RECREATE");
  fout->WriteTObject(hc[1],"analyzed_1","TObject::kOverwrite");
  fout->Close();









}
