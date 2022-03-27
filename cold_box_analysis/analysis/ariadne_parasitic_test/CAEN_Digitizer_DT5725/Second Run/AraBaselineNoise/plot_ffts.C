void setHistograms(TH1D *h){

  h->GetYaxis()->SetTitle("Magnitude");
  h->GetXaxis()->SetTitle("Frequency (MHz)");
  h->GetXaxis()->SetRangeUser(0,125);

}

void plot_ffts(){

  TFile *f1 = new TFile("graphs/ffts_no_filter.root","READ");


  THStack *hs1 = (THStack*)f1->Get("hstack1");
  THStack *hs2 = (THStack*)f1->Get("hstack2");

  TList *ml1 = (TList*)hs1->GetHists();
  TList *ml2 = (TList*)hs2->GetHists();


  TH1D *h1t = (TH1D*)ml1->At(0);
  TH1D *h1c = (TH1D*)ml1->At(1);
  TH1D *h2t = (TH1D*)ml2->At(0);
  TH1D *h2c = (TH1D*)ml2->At(1);

  setHistograms(h1t);
  setHistograms(h1c);
  setHistograms(h2t);
  setHistograms(h2c);
  
  TCanvas *c = new TCanvas("c","c",1920,0,1920,1080);
  h1c->Draw();
  h2c->Draw("SAME");
}
