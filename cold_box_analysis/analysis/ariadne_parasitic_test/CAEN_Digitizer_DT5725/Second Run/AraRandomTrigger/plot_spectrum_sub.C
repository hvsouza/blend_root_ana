void setHistograms(TH1D *h){

  h->GetYaxis()->SetTitle("Magnitude");
  h->GetXaxis()->SetTitle("Frequency (MHz)");
  h->GetXaxis()->SetRangeUser(0,90);

}

void plot_spectrum_sub(){


  TFile *fnoise = new TFile("graphs/ffts_no_filter.root","READ");
  TFile *fsignal = new TFile("../AraBaselineNoise/graphs/ffts_no_filter.root","READ");

  THStack *hnoise1 = (THStack*)fnoise->Get("hstack1");
  THStack *hnoise2 = (THStack*)fnoise->Get("hstack2");

  THStack *hsignal1 = (THStack*)fsignal->Get("hstack1");
  THStack *hsignal2 = (THStack*)fsignal->Get("hstack2");

  TList *mln1 = (TList*)hnoise1->GetHists();
  TList *mln2 = (TList*)hnoise2->GetHists();

  TList *mls1 = (TList*)hsignal1->GetHists();
  TList *mls2 = (TList*)hsignal2->GetHists();


  TH1D *hn1c = (TH1D*)mln1->At(1);
  TH1D *hn2c = (TH1D*)mln2->At(1);

  TH1D *hs1c = (TH1D*)mls1->At(1);
  TH1D *hs2c = (TH1D*)mls2->At(1);
  setHistograms(hn1c);
  setHistograms(hn2c);
  setHistograms(hs1c);
  setHistograms(hs2c);

  hn1c->SetLineColor(kBlack);
  hn2c->SetLineColor(kBlack);

  hn1c->SetTitle("C1XArapuca noise (SiPMs off)");
  hn2c->SetTitle("C2XArapuca noise (SiPMs off)");
  hs1c->SetTitle("C1XArapuca signal");
  hs2c->SetTitle("C2XArapuca signal");
  
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  gPad->SetLogy(1);
  hs1c->Draw();
  hn1c->Draw("SAME");
  hs1c->GetYaxis()->SetRangeUser(4e6,3e8);
  
  TCanvas *c2 = new TCanvas("c2","c2",1920,0,1920,1080);
  gPad->SetLogy(1);
  hs2c->Draw();
  hn2c->Draw("SAME");
  hs2c->GetYaxis()->SetRangeUser(4e6,3e8);

  TH1D *hsub1 = (TH1D*)hs1c->Clone("hsub1");
  TH1D *hsub2 = (TH1D*)hs2c->Clone("hsub2");

  hsub1->Add(hs1c,hn1c,1,-1);
  hsub2->Add(hs2c,hn2c,1,-1);

  setHistograms(hsub1);
  setHistograms(hsub2);
  hsub1->SetTitle("C1XArapuca subtracted");
  hsub2->SetTitle("C2XArapuca subtracted");
  
  TCanvas *csub1 = new TCanvas("csub1","csub1",1920,0,1920,1080);
  gPad->SetLogy(1);
  hsub1->Draw();
  hsub1->GetYaxis()->SetRangeUser(1e2,7e8);
  
  TCanvas *csub2 = new TCanvas("csub2","csub2",1920,0,1920,1080);
  gPad->SetLogy(1);
  hsub2->Draw();
  hsub2->GetYaxis()->SetRangeUser(1e2,7e8);
  c1->BuildLegend();
  c2->BuildLegend();
  csub1->BuildLegend();
  csub2->BuildLegend();
}

