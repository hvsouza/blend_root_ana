#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_wvfs_v5_ch2_large(){
  // gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch5";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0"){minInt = 2600; maxInt = 3500;}
  if(mychannel == "Ch2"){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch3"){minInt = 2600; maxInt = 3500;}
  if(mychannel == "Ch4"){minInt = 2600; maxInt = 3500;}
  if(mychannel == "Ch5"){minInt = 2600; maxInt = 3500;}
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

  vector<string> files = {"./run9_cosmic_all_trigger_xArapuca_v5_Ch2_400ADCs",
                          "./run10_cosmic_all_trigger_xArapuca_v5_Ch2_2000ADCs"};

  vector<Double_t> volts = {1,2};
  Double_t saturation_level;

  Int_t n = volts.size();
  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);

  TH2D *h2p = new TH2D("h2p","h2p",4000,-5,2.3e6,2000,-100,14000);
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  c1->DrawFrame(0,-.1,20000,1.1);

  Double_t sphe = 1;
  vector<vector<Double_t>> group_amp = {{500,600}, {800,900}, {1000,1100}, {1500,1600}, {1800,1900}, {2500,2600}, {3200,3300}, {4000,4200}, {7200,7300}, {8000,8200}, {10000,10200}, {11000,12000}};
  Int_t ngroups = group_amp.size();
  vector<TH1D*> h(ngroups);
  vector<Int_t> nev(ngroups);
  THStack *hs = new THStack("hs", "hs");

  for(Int_t k = 0; k < ngroups; k++){
    h[k] = new TH1D(Form("h%d",k),Form("h%d",k),2500, 0, 10000);
  }
  ANALYZER *a = nullptr;
  for(Int_t i = 0; i<n; i++){
    cout << "File: " << files[i] << endl;
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z%.2f",volts[i]));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    a = z[i];
    Int_t kch = z[i]->kch;
    saturation_level = saturations[z[i]->getIdx()];
    sphe = sphes[z[i]->getIdx()];
    for(Int_t j = 0; j < a->nentries; j++){
      a->getWaveform(j);
      a->applyDenoise(16);
      a->integrate(minInt,maxInt,0.1,false);
      Double_t max = a->temp_max;
      h2p->Fill(a->temp_charge,a->temp_max);
      for(Int_t k = 0; k < ngroups; k++){
        if(max >= group_amp[k][0] && max <= group_amp[k][1]){
          for(Int_t l = 0; l < memorydepth; l++){
            h[k]->AddBinContent(l+1,a->ch[kch].wvf[l]);
          }
          nev[k] += 1;
        }
      }
    }

    h2p->Draw("colz");
  }

  TCanvas *c2 = new TCanvas("c2", "c2",1920,0,1920,1080);
  c2->DrawFrame(0,-0.1,10000,1.1);

  c2->cd();

  for(Int_t k = 0; k < ngroups; k++){
    h[k]->Scale(1./nev[k]);
    Double_t max = h[k]->GetMaximum();
    h[k]->Scale(1./max);
    stringstream stream;
    stream << std::fixed << std::setprecision(1) << max;
    string correct_max = stream.str();
    string hist_name = Form("%s - Max = %s (%d entries)",devices[a->getIdx()].c_str(), correct_max.c_str(), nev[k]);
    h[k]->SetTitle(hist_name.c_str());
    h[k]->SetLineWidth(2);
    hs->Add(h[k],"hist");
  }
  hs->GetXaxis()->SetTitle("Time (ns)");
  hs->GetYaxis()->SetTitle("Amplitude (A.U.)");
  hs->Draw("nostack plc same");
  c2->BuildLegend();


}
#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_wvfs_Argon2x2_large(){
  // gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch0";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0"){minInt = 2600; maxInt = 3500;}
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

  vector<string> files = {"./run14_cosmic_all_trigger_Argon2x2_2000ADCs"};
  vector<Double_t> volts = {1};
  Double_t saturation_level;

  Int_t n = volts.size();
  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);

  TH2D *h2p = new TH2D("h2p","h2p",4000,-5,2.3e6,2000,-100,14000);
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  c1->DrawFrame(0,-.1,20000,1.1);

  Double_t sphe = 1;
  vector<vector<Double_t>> group_amp = {{2500,2600}, {4000,4200}, {7000,7200}, {10000,12000}};
  Int_t ngroups = group_amp.size();
  vector<TH1D*> h(ngroups);
  vector<Int_t> nev(ngroups);
  THStack *hs = new THStack("hs", "hs");

  for(Int_t k = 0; k < ngroups; k++){
    h[k] = new TH1D(Form("h%d",k),Form("h%d",k),2500, 0, 10000);
  }
  ANALYZER *a = nullptr;
  for(Int_t i = 0; i<n; i++){
    cout << "File: " << files[i] << endl;
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z%.2f",volts[i]));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    a = z[i];
    Int_t kch = z[i]->kch;
    saturation_level = saturations[z[i]->getIdx()];
    sphe = sphes[z[i]->getIdx()];
    for(Int_t j = 0; j < a->nentries; j++){
      a->getWaveform(j);
      a->applyDenoise(16);
      a->integrate(minInt,maxInt,0.1,false);
      Double_t max = a->temp_max;
      h2p->Fill(a->temp_charge,a->temp_max);
      for(Int_t k = 0; k < ngroups; k++){
        if(max >= group_amp[k][0] && max <= group_amp[k][1]){
          for(Int_t l = 0; l < memorydepth; l++){
            h[k]->AddBinContent(l+1,a->ch[kch].wvf[l]);
          }
          nev[k] += 1;
        }
      }
    }

    h2p->Draw("colz");
  }

  TCanvas *c2 = new TCanvas("c2", "c2",1920,0,1920,1080);
  c2->DrawFrame(0,-0.1,10000,1.1);

  c2->cd();

  for(Int_t k = 0; k < ngroups; k++){
    h[k]->Scale(1./nev[k]);
    Double_t max = h[k]->GetMaximum();
    h[k]->Scale(1./max);
    stringstream stream;
    stream << std::fixed << std::setprecision(1) << max;
    string correct_max = stream.str();
    string hist_name = Form("%s - Max = %s (%d entries)",devices[a->getIdx()].c_str(), correct_max.c_str(), nev[k]);
    h[k]->SetTitle(hist_name.c_str());
    h[k]->SetLineWidth(2);
    hs->Add(h[k],"hist");
  }
  hs->GetXaxis()->SetTitle("Time (ns)");
  hs->GetYaxis()->SetTitle("Amplitude (A.U.)");
  hs->Draw("nostack plc same");
  c2->BuildLegend();


}
