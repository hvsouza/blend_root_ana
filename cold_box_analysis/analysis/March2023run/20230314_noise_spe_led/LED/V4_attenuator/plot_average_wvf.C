#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_average_wvf(){
  gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch4.";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0."){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch2."){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch3."){minInt = 10300; maxInt = 12000;}
  if(mychannel == "Ch4."){minInt = 10300; maxInt = 10720;}
  vector<string> devices = {"miniArapuca 37V (A1ch1)","none", "none", "xArapuca v4 ch2 (DCem1.2)", "xArapuca v4 ch1 (DCem1.2)", "xArapuca v5 ch2 (DCemArgon4)", "xArapuca v5 ch1 (DCemArgon4)"};

  vector<Double_t> saturations = {12000, 12000, 12600, 12600, 12000, 12000, 12000, 100};
  // vector<Double_t> sphes = {7100,400, 7100, 5800*9.97/515.25, 5800, 0, 0, 0};
  vector<Double_t> sphes = {7100,400, 7100, 1, 5800, 0, 0, 0};

  vector<string> files = {"run32_v4_5dB_ch4_10dB_ch5_365nm_20ns_3V4", "run33_v4_5dB_ch4_10dB_ch5_365nm_20ns_3V6", "run34_v4_5dB_ch4_10dB_ch5_365nm_20ns_3V8", "run35_v4_5dB_ch4_10dB_ch5_365nm_20ns_4V0", "run36_v4_5dB_ch4_10dB_ch5_365nm_20ns_4V2", "run37_v4_5dB_ch4_10dB_ch5_365nm_20ns_4V4", "run38_v4_5dB_ch4_10dB_ch5_365nm_20ns_4V6", "run39_v4_5dB_ch4_10dB_ch5_365nm_20ns_4V8", "run40_v4_5dB_ch4_10dB_ch5_365nm_20ns_5V0", "run41_v4_5dB_ch4_10dB_ch5_365nm_20ns_5V2", "run42_v4_5dB_ch4_10dB_ch5_365nm_20ns_5V4", "run43_v4_5dB_ch4_10dB_ch5_365nm_20ns_5V6", "run44_v4_5dB_ch4_10dB_ch5_365nm_20ns_5V8", "run45_v4_5dB_ch4_10dB_ch5_365nm_20ns_6V0", "run46_v4_5dB_ch4_10dB_ch5_365nm_20ns_6V5", "run47_v4_5dB_ch4_10dB_ch5_365nm_20ns_7V0"};
  vector<Double_t> volts = {3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.5, 7.0};
  Double_t saturation_level;

  Int_t n = volts.size();
  vector<Int_t> range_files(2);
  if(mychannel == "Ch0.") range_files = {2,24};
  if(mychannel == "Ch3.") range_files = {1, 15};
  if(mychannel == "Ch4.") range_files = {0,15};
  if(mychannel == "Ch5.") range_files = {1,12};



  Int_t nbins;
  Double_t dtime = 4;
string conc = "/analyzed.root"; vector<ANALYZER*> z(n);

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
    saturation_level = saturations[z[i]->getIdx()];
    sphe = sphes[z[i]->getIdx()];

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
