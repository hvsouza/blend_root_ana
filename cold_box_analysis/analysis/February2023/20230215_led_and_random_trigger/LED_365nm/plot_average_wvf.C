#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_average_wvf(){
  // gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch2";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0"){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch2"){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch3"){minInt = 10300; maxInt = 12000;}
  if(mychannel == "Ch4"){minInt = 10300; maxInt = 10720;}
  vector<string> devices = {
                            "Argon2x2",
                            "Argon4",
                            "DCemV1.0",
                            "xArapuca V4 Ch2",
                            "xArapuca V4 Ch1",
                            "xArapuca V5 Ch2",
                            "xArapuca V5 Ch1",
                            "DCemVD1.2",};

  vector<Double_t> saturations = {12000, 12000, 12600, 12600, 12000, 12000, 12000, 100};
  // vector<Double_t> sphes = {7100,400, 7100, 5800*9.97/515.25, 5800, 0, 0, 0};
  vector<Double_t> sphes = {7100,400, 7100, 1, 5800, 0, 0, 0};

  vector<string> files = {"run44_all_devices_365nm_20ns_3V2", "run45_all_devices_365nm_20ns_3V4", "run46_all_devices_365nm_20ns_3V6", "run47_all_devices_365nm_20ns_3V8", "run48_all_devices_365nm_20ns_4V0", "run49_all_devices_365nm_20ns_4V2", "run50_all_devices_365nm_20ns_4V4", "run51_all_devices_365nm_20ns_4V6", "run52_all_devices_365nm_20ns_4V8", "run53_all_devices_365nm_20ns_5V0", "run54_all_devices_365nm_20ns_5V2", "run55_all_devices_365nm_20ns_5V4", "run56_all_devices_365nm_20ns_5V6", "run57_all_devices_365nm_20ns_5V8", "run58_all_devices_365nm_20ns_6V0", "run59_all_devices_365nm_20ns_6V5", "run60_all_devices_365nm_20ns_7V0", "run61_all_devices_365nm_20ns_7V5", "run62_all_devices_365nm_20ns_8V0", "run63_all_devices_365nm_20ns_8V5", "run64_all_devices_365nm_20ns_9V0", "run65_all_devices_365nm_20ns_10V0", "run66_all_devices_365nm_20ns_11V0", "run67_all_devices_365nm_20ns_12V5", "run68_all_devices_365nm_20ns_15V0", "run69_all_devices_365nm_20ns_17V5", "run70_all_devices_365nm_20ns_20V0", "run71_all_devices_365nm_20ns_22V5", "run72_all_devices_365nm_20ns_25V0", "run73_all_devices_365nm_20ns_30V0"};
  vector<Double_t> volts = {3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.0, 11.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 30.0};
  Double_t saturation_level;

  Int_t n = volts.size();
  vector<Int_t> range_files(2);
  if(mychannel == "Ch0") range_files = {2,24};
  if(mychannel == "Ch2") range_files = {1,15};
  if(mychannel == "Ch3") range_files = {0,15};
  if(mychannel == "Ch4") range_files = {1, 8};
  if(mychannel == "Ch5") range_files = {1,12};

  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);

  vector<TH1D*> h(n);
  THStack *hs = new THStack("hs", "hs");
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  c1->DrawFrame(0,-.1,20000,1.1);

  Double_t sphe = 1;

  ANALYZER a;
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
    hs->GetXaxis()->SetTitle("Time (ns)");
    hs->GetYaxis()->SetTitle("Amplitude (A.U.)");
  }

  hs->Draw("nostack plc same");
  c1->BuildLegend();

}
