void DCR_analysis_V1(){

  Int_t nfiles = 6;
  vector<TFile*> f1(nfiles);
  vector<string> fname = {"run2_switching_on_argon2x2_argon4_DCem1dot2VD_47V_DCem1dot0VD_36V",
                          "run3_switching_on_argon2x2_argon4_DCem1dot2VD_47V_DCem1dot0VD_36V_with_SSP",
                          "run4_switching_on_argon2x2_argon4_DCem1dot0VD_36V_with_SSP",
                          "run5_switching_on_argon2x2_argon4_DCem1dot2VD_47VD_DCem1dot0VD_36V_with_SSP__xa_V4",
                          "run6_switching_on_argon2x2_argon4_DCem1dot2VD_47VD_DCem1dot0VD_36V_with_SSP__xa_V4_and_xa_V5",
                          "run7_switching_on_argon2x2_argon4_DCem1dot2VD_47VD_DCem1dot0VD_36V_with_SSP_xa_V5"
                          };

  vector<string> description = {"miniArapucas + flex + V1 (VD)",
                                "miniArapucas + flex + V1 (VD-SSP)",
                                "miniArapucas +  V1 (VD-SSP) ",
                                "miniArapucas + flex + V1 (VD) + v4",
                                "miniArapucas + flex + V1 (VD) + v4 +v5",
                                "miniArapucas + flex + V1 (VD) + v5",
                                };

  vector<TH1D*> h(nfiles);
  THStack *hs = new THStack("hs","hs");
  Double_t total_entries = 0;
  Double_t total_time = 500*1e-3;
  Double_t total_sipms_area = 36*80;
  Double_t dcr = 0;

  for(Int_t i = 0; i < nfiles; i++){
    fname[i] = fname[i] + "/sphe_histograms_Ch2.root";
    f1[i] = new TFile(fname[i].c_str(),"READ");
    h[i] = (TH1D*)f1[i]->Get("analyzed_2");
    h[i]->Rebin(200);
    Double_t minval = h[i]->GetBinCenter(1);
    Double_t binwidth = h[i]->GetBinWidth(1);
    Int_t binCut = round(500 - minval)/binwidth + 1;
    total_entries = h[i]->Integral(binCut,h[i]->GetNbinsX());
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
