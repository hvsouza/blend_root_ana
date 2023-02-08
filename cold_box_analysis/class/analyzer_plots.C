// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer_plots.C
// Created: 2023-02-08
// ________________________________________ //
//

#include "MYCODES.h"

/**
 *
 * This class holds the methods for ploting for analyzer */
void ANALYZER::sample_plot(Int_t myevent = 0, string opt = "", Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
  if (opt == "") opt = plot_opt;
  bool state = getWaveformHard(myevent,factor);
  if (!state) return;
  applyMovingAverage(mafilter);
  applyDenoise(filter);

  drawGraph(opt,n_points,&time[0],&ch[kch].wvf[0]);
}
void showWaveform(Int_t maxevent = 0, Int_t filter = 0, Int_t dt = 0){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TCanvas *c1 = new TCanvas("c1");

  for(Int_t i = 0; i < maxevent; i++){
    sample_plot(i,"AL",filter);
    printf("\rEvent %d", i);
    fflush(stdout);
    c1->Modified();
    c1->Update();
    if (gSystem->ProcessEvents())
      break;
    if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
  }
}
void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="", Double_t factor = 1){

  Int_t nbinsx = (xmax-xmin)/dtime;
  if(!hpers) hpers = new TH2D("hpers","hpers",nbinsx,xmin,xmax,nbins,ymin,ymax);
  else{
    hpers->Reset();
    hpers->SetBins(nbinsx, xmin, xmax, nbins, ymin, ymax);
  }
  if(!cpers) cpers = new TCanvas(Form("cpers_%s", myname.c_str()), Form("cpers_%s", myname.c_str()),1920,0,1920,1080);
  else{cpers->cd();}



  getSelection(cut);
  Int_t nev = lev->GetN();
  Int_t iev = 0;
  for(Int_t i = 0; i < nev; i++){
    iev = lev->GetEntry(i);
    getWaveform(iev,kch);
    applyDenoise(filter);

    for (int j = 0; j < n_points; j++) {
      hpers->Fill(j*dtime,ch[kch].wvf[j]*factor);
    }

  }

  hpers->Draw("colz");
  hpers->GetXaxis()->SetTitle("Time (ns)");
  hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

}


void drawGraph(string opt = "", Int_t n = memorydepth, Double_t* x = nullptr, Double_t* y = nullptr){
  if (opt == "") opt = plot_opt;
  if (x == nullptr) x = time;
  if (y == nullptr) y = ch[kch].wvf;
  gwvf = new TGraph(n,x,y);
  gwvf->Draw(opt.c_str());
  gwvf->GetXaxis()->SetTitle(xlabel.c_str());
  gwvf->GetYaxis()->SetTitle(ylabel.c_str());

  if(ymax!=0 && ymin!=0){
    gwvf->GetYaxis()->SetRangeUser(ymin,ymax);
  }
  gwvf->GetXaxis()->SetRangeUser(xmin,xmax);
  gwvf->SetEditable(kFALSE);
}

void ANALYZER::averageFFT(Int_t maxevent = 0, string selection = "", bool inDecibel = false, Double_t filter = 0){
  if (maxevent==0) {
    maxevent = nentries;
  }
  getSelection(selection);
  Int_t nev = lev->GetN();
  if (maxevent < nev) {
    nev = maxevent;
  }
  Int_t iev = 0;
  hfft[kch] = (TH1D*)w->hfft->Clone(Form("hfft_%s_ch%d",myname.c_str(),kch));
  hfft[kch]->Reset();
  Int_t total = 0;
  for(Int_t i = 0; i < nev; i++){
    iev = lev->GetEntry(i);
    getWaveform(iev,kch);
    applyDenoise(filter);
    getFFT();
    for (Int_t j = 0; j < memorydepth/2; j++) hfft[kch]->AddBinContent(j+1,w->hfft->GetBinContent(j+1));
    total++;
  }
  hfft[kch]->Scale(1./total);
  hfft[kch]->SetEntries(total);

  if (inDecibel){
    w->convertDecibel(hfft[kch]);
    hfft[kch]->GetYaxis()->SetTitle("Magnitude (dB)");
  }
}
void ANALYZER::showFFT(Int_t naverage = 10, Int_t maxevent = 0, Int_t dt = 0, bool inDecibel = false){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TH1D *h = (TH1D*)w->hfft->Clone("h");
  vector<TH1D *> hfft(naverage);
  Int_t k = 0;
  TCanvas *c1 = new TCanvas("c1");
  for(Int_t i = 0; i < maxevent; i++){
    getWaveform(i,kch);
    getFFT();
    if (inDecibel) w->convertDecibel();
    if (i < naverage) {
      hfft[k] = (TH1D*)w->hfft->Clone(Form("h%d",k));
      for (Int_t j = 0; j < memorydepth/2; j++) h->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
      if (naverage == 1) h->Draw();
    }
    else{
      if (k < naverage) {
        for (Int_t j = 0; j < memorydepth/2; j++) {
          h->AddBinContent(j+1,-hfft[k]->GetBinContent(j+1));
          hfft[k]->SetBinContent(j+1,w->hfft->GetBinContent(j+1));
          h->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
        }
        printf("\rEvent %d", i);
        fflush(stdout);
        c1->Modified();
        c1->Update();
        if (gSystem->ProcessEvents())
          break;
        if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
      }
      else{
        if (i == naverage){
          h->Draw();
        }
        k = 0;
      }
    }

    k++;
    if(naverage == 1) k = 0;
  }
  void ANALYZER::debugSPE(Int_t event, Int_t moving_average, Int_t n_moving, Int_t shift, vector<Double_t> xrange, vector<Double_t> yrange){


  }

}
