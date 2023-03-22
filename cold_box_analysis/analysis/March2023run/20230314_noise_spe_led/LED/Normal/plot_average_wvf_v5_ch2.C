#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_average_wvf_v5_ch2(){
  gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch5.";
  vector<string> devices = {"miniArapuca 37V (A1ch1)","none", "none", "xArapuca v4 ch2 (DCem1.2)", "xArapuca v4 ch1 (DCem1.2)", "xArapuca v5 ch2 (DCemArgon4)", "xArapuca v5 ch1 (DCemArgon4)"};

  vector<Double_t> saturations = {12000, 12000, 12600, 12600, 12000, 12000, 12000, 100};

  vector<string> files = {"run7_argon2x2_v4_v5_365nm_20ns_3V2", "run8_argon2x2_v4_v5_365nm_20ns_3V4", "run9_argon2x2_v4_v5_365nm_20ns_3V6", "run10_argon2x2_v4_v5_365nm_20ns_3V8", "run11_argon2x2_v4_v5_365nm_20ns_4V0", "run12_argon2x2_v4_v5_365nm_20ns_4V2", "run13_argon2x2_v4_v5_365nm_20ns_4V4", "run14_argon2x2_v4_v5_365nm_20ns_4V6", "run15_argon2x2_v4_v5_365nm_20ns_4V8", "run16_argon2x2_v4_v5_365nm_20ns_5V0", "run17_argon2x2_v4_v5_365nm_20ns_5V2", "run18_argon2x2_v4_v5_365nm_20ns_5V4", "run19_argon2x2_v4_v5_365nm_20ns_5V6", "run20_argon2x2_v4_v5_365nm_20ns_5V8", "run21_argon2x2_v4_v5_365nm_20ns_6V0", "run22_argon2x2_v4_v5_365nm_20ns_6V5", "run23_argon2x2_v4_v5_365nm_20ns_7V0", "run24_argon2x2_v4_v5_365nm_20ns_7V5", "run25_argon2x2_v4_v5_365nm_20ns_8V0", "run26_argon2x2_v4_v5_365nm_20ns_8V5", "run27_argon2x2_v4_v5_365nm_20ns_9V0", "run28_argon2x2_v4_v5_365nm_20ns_9V5", "run29_argon2x2_v4_v5_365nm_20ns_10V0", "run30_argon2x2_v4_v5_365nm_20ns_11V0", "run31_argon2x2_v4_v5_365nm_20ns_12V0"};

  vector<Double_t> volts = {3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0};

  Int_t n = volts.size();
  vector<Int_t> range_files(2);
  if(mychannel == "Ch0.") range_files = {2,24};
  if(mychannel == "Ch3.") range_files = {1, 15};
  if(mychannel == "Ch4.") range_files = {1,n};
  if(mychannel == "Ch5.") range_files = {1,n-5};
  if(mychannel == "Ch6.") range_files = {1,11};



  Int_t nbins;
  Double_t dtime = 4;
  string conc = "/analyzed.root";
  vector<ANALYZER*> z(n);

  vector<TH1D*> h(n);
  THStack *hs = new THStack("hs", "hs");
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  c1->DrawFrame(0,-.1,20000,1.1);

  Double_t sphe = 1;
  TLegend *lg = new TLegend(0.8,0.8,1,1);

  for(Int_t i = range_files[0]; i<range_files[1]; i++){
    cout << "File: " << files[i] << endl;
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z%.2f",volts[i]));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    Int_t kch = z[i]->kch;

    z[i]->averageWaveform(500);
    h[i] = (TH1D*)z[i]->h->Clone(Form("h%d",i));
    Double_t max = h[i]->GetMaximum();
    h[i]->Scale(1./max);
    stringstream stream;
    stream << std::fixed << std::setprecision(1) << max;
    string correct_max = stream.str();
    string hist_name = Form("%s - LED %.1f V - Max = %s",devices[z[i]->getIdx()].c_str(),volts[i], correct_max.c_str());
    h[i]->SetTitle(hist_name.c_str());
    h[i]->SetLineWidth(2);
    hs->Add(h[i],"hist");
    lg->AddEntry(h[i], hist_name.c_str(),"l");
    hs->GetXaxis()->SetTitle("Time (ns)");
    hs->GetYaxis()->SetTitle("Amplitude (A.U.)");
  }

  hs->Draw("nostack plc same");
  lg->Draw();

}
