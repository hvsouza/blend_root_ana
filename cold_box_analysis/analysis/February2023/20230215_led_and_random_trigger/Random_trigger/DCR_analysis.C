void DCR_analysis(){

  Int_t nfiles = 5;
  vector<TFile*> f1(nfiles);

  vector<string> fname = {
                          "run74_all_devices",
                          "run75_argon2x2_argon4_DCemVD1dot0_DCemVD1dot2_xArapucaV5",
                          "run76_argon2x2_argon4_DCemVD1dot0_DCemVD1dot2",
                          "run77_argon2x2_argon4_DCemVD1dot0",
                          "run78_argon2x2_argon4"
                          };

  vector<string> description = {
                                "miniArapucas + flex + V1 (VD-SSP) + v4 +v5",
                                "miniArapucas + flex + V1 (VD-SSP) + v5",
                                "miniArapucas + flex + V1 (VD-SSP)",
                                "miniArapucas +  V1 (VD-SSP) ",
                                "miniArapucas"
                                };

  vector<TH1D*> h(nfiles);
  THStack *hs = new THStack("hs","hs");
  Double_t total_entries = 0;
  Double_t total_time = 500*1e-3;
  Double_t total_sipms_area = 36*20;
  Double_t dcr = 0;

  for(Int_t i = nfiles-1; i >= 0; i--){
    fname[i] = fname[i] + "/sphe_histograms_Ch0.root";
    f1[i] = new TFile(fname[i].c_str(),"READ");
    h[i] = (TH1D*)f1[i]->Get("analyzed_0");
    h[i]->Rebin(200);
    total_entries = h[i]->GetEntries();
    dcr = total_entries/total_time/total_sipms_area;
    stringstream stream;
    stream << std::fixed << std::setprecision(2) << dcr;
    description[i] = description[i] + ", DCR = " + stream.str() + " Hz/mm^{2}";
    h[i]->SetTitle(description[i].c_str());
    h[i]->SetLineWidth(3);
    // h[i]->Scale(1./h[i]->GetBinWidth(1));
    h[i]->Scale(1.,"width");
    hs->Add(h[i], "hist");
  }

  TCanvas *c = new TCanvas("c", "c",1920,0,1920,1080);

  hs->Draw("nostack plc");
  hs->GetXaxis()->SetTitle("Charge (ADC*ns)");
  hs->GetYaxis()->SetTitle("counts/bin width");
  c->BuildLegend();
}
