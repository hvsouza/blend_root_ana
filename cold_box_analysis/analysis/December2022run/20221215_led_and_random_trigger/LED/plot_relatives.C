// ________________________________________ //
// Author: Henrique Souza
// Filename: plot_relatives.C
// Created: 2023-01-31
// ________________________________________ //
//



void plot_relatives(){
  // gStyle->SetPalette(kViridis);
  vector<Int_t> th = {200,300,400,600,800,1000};
  Int_t n = th.size();

  vector<TFile *> f(n);
  vector<TCanvas *> c(n);
  vector<TH1D *> h(n);

  vector<TFile *> fa(n);
  vector<TCanvas *> ca(n);
  vector<TH1D *> ha(n);
  THStack *hs = new THStack("hs","hs");
  THStack *hsa = new THStack("hsa","hsa");

  for (Int_t i = 0; i < n; i++) {
    f[i] = new TFile(Form("graphs/TOT_relative_th%d.root",th[i]), "READ");
    fa[i] = new TFile(Form("graphs/TOT_relative_amp_th%d.root",th[i]), "READ");
    c[i] = (TCanvas*)f[i]->Get("chdev");
    ca[i] = (TCanvas*)fa[i]->Get("chdev");
    h[i] = (TH1D*)c[i]->GetPrimitive("hdev");
    h[i]->SetTitle(Form("Threshold: %d", th[i]));
    h[i]->SetLineWidth(3);
    hs->Add(h[i]);

    ha[i] = (TH1D*)ca[i]->GetPrimitive("hdev");
    ha[i]->SetTitle(Form("Threshold: %d", th[i]));
    ha[i]->SetLineWidth(3);
    hsa->Add(ha[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  hs->Draw("nostack PLC");
  c1->BuildLegend();

  TCanvas *c2 = new TCanvas("c2","c2");
  hsa->Draw("nostack PLC");
  c2->BuildLegend();
}
