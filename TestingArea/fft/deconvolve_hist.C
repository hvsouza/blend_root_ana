#define memorydepth 400
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void deconvolve_hist(){
  TFile *f = new TFile("sphe_histograms_darkCount_Ch5.root","READ");
  TH1D *h = (TH1D*)f->Get("analyzed_5");
  h->Rebin(125);
  Int_t nbins = h->GetNbinsX();

  WIENER w("w",nbins);
  WIENER g("g",nbins);

  w.setFromHist(h);
  g.setFromHist(h);


  w.fft(h);
  g.fillPlane();

  TH1D *hg = (TH1D*)h->Clone("hg");
  hg->Reset();
  TF1 *fga = new TF1("fga","TMath::Gaus(x,[0],[1])",-20000,80000);
  fga->SetParameters(0,600);
  for (Int_t i = 0; i < nbins; i++) {
    hg->SetBinContent(i+1,fga->Eval(i*hg->GetBinWidth(1)));
  }
  g.setFromHist(hg);
  g.fft(g.hwvf);

  // g.setFilter(0.05, "gaus");
  // g.apply_filter();
  // g.backfft(g);


  TCanvas *c1 = new TCanvas("c1");
  c1->Divide(2,1);
  c1->cd(1);
  TH1D *h1 = (TH1D*)w.hwvf->Clone("h1");
  TH1D *h2 = (TH1D*)w.hfft->Clone("h2");
  h1->Draw();
  h1->GetYaxis()->SetRangeUser(0,h1->GetMaximum());
  c1->cd(2);
  h2->Draw();

  TCanvas *c2 = new TCanvas("c2");
  c2->Divide(2,1);
  c2->cd(1);
  g.hwvf->Draw();
  c2->cd(2);
  g.hfft->Draw();


  w.deconvolve(w, g,0.046);

  TCanvas *c3 = new TCanvas("c3");
  c3->Divide(2,1);
  c3->cd(1);
  w.hwvf->Draw();
  w.hwvf->GetYaxis()->SetRangeUser(0,w.hwvf->GetMaximum());
  h1->SetLineColor(kRed);
  // h1->Draw("SAME");
  c3->cd(2);
  w.hfft->Draw();


  TCanvas *c5 = new TCanvas();
  TSpectrum *s = new TSpectrum();
  Double_t source[nbins];
  Double_t response[nbins];
  for (Int_t i = 0; i < nbins; i++) {
    source[i] = h->GetBinContent(i+1);
    response[i] = g.hwvf->GetBinContent(i+1);
  }
  s->Deconvolution(source,response,nbins,1,1,1);
  TH1D *hs = (TH1D*)h->Clone("hs");
  hs->Reset();
  for (Int_t i = 0; i < nbins; i++) hs->SetBinContent(i + 1,source[i]);
  hs->SetLineColor(kRed);
  // h->Draw();
  hs->Draw();




}
