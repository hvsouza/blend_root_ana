void DCR_analysis(){

  vector<string> fname = {
                          "./run0_argon2x2_v4_v5_starting_of_day",
                          "./run1_argon2x2_starting_of_day",
                          "./run19_argon2x2_v4_v5_cathode_off",
                          "./run20_argon2x2_cathode_off",
                          "./run21_argon2x2_v4_v5_cathode_off"
                          };

  vector<string> description = {
                                "V4 and V5 on (start of the day)",
                                "No PoF",
                                "V4 and V5 (HV ON)",
                                "PoF off (HV OFF)",
                                "V4 and V5 on (HV OFF)"
                                };
  Int_t nfiles = fname.size();
  vector<TFile*> f1(nfiles);
  vector<TH1D*> h(nfiles);
  THStack *hs = new THStack("hs","hs");
  Double_t total_entries = 0;
  Double_t total_time = 500*1e-3;
  Double_t total_sipms_area = 36*20;
  Double_t dcr = 0;

  for(Int_t i = 0; i < nfiles; i++){
    fname[i] = fname[i] + "/sphe_histograms_Ch0.root";
    f1[i] = new TFile(fname[i].c_str(),"READ");
    h[i] = (TH1D*)f1[i]->Get("analyzed");
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
