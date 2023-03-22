#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_wvfs_v4_ch1(){
  // gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch4.";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0."){minInt = 2600; maxInt = 3500;}
  if(mychannel == "Ch2."){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch3."){minInt = 2600; maxInt = 3500;}
  if(mychannel == "Ch4."){minInt = 2600; maxInt = 3500;}
  vector<string> devices = {
                            "Argon2x2",
                            "Argon4",
                            "DCemV1.0",
                            "xArapuca V4 Ch2",
                            "xArapuca V4 Ch1",
                            "xArapuca V5 Ch2",
                            "xArapuca V5 Ch1",
                            "DCemVD1.2",};

  vector<Double_t> saturations = {12000, 12000, 12600, 12600, 13400, 12000, 13000, 100};
  vector<Double_t> sphes = {7100,400, 7100, 1, 4200, 0, 2500};
  vector<Double_t> sphes_amp = {44.7, 0, 0, 18.2, 26.9, 4.2, 29.6};

  vector<string> files = {
                          "run24_argon2x2_v4_v5_trigger_Ch4_700ADC",
                          "run25_argon2x2_v4_v5_trigger_Ch4_2500ADC"
                          };

  Int_t n = files.size();
  Double_t saturation_level;

  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);

  TH2D *h2p = new TH2D("h2p","h2p",500,0,5000,200,0,1);
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  c1->DrawFrame(0,-.1,20000,1.1);

  Double_t spe = 1;
  Double_t spe_amp= 1;
  vector<vector<Double_t>> group_amp = {
                                       {50,100},
                                       {100,200},
                                       {300,400},
                                       {400,600},
                                       };
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
    z[i] = new ANALYZER(Form("z%d",i));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    a = z[i];
    Int_t kch = z[i]->kch;
    saturation_level = saturations[z[i]->getIdx()];
    spe = sphes[z[i]->getIdx()];
    spe_amp = sphes_amp[z[i]->getIdx()];
    for(Int_t j = 0; j < a->nentries; j++){
      a->getWaveform(j);
      // a->applyDenoise(16);
      a->integrate(minInt,maxInt,0.1,false);
      Double_t max = a->temp_max;
      Double_t charge = a->temp_charge;
      if(max > saturation_level) continue;
      for(Int_t k = 0; k < ngroups; k++){
        if(max/spe_amp >= group_amp[k][0] && max/spe_amp <= group_amp[k][1]){
          for(Int_t l = 0; l < a->n_points; l++){
            h[k]->AddBinContent(l+1,a->ch[kch]->wvf[l]);
          }
          nev[k] += 1;
        }

      }
      a->integrate(2800,3050);
      Double_t cfast = a->temp_charge;
      if(max > 700) h2p->Fill(charge/spe,cfast/charge);
    }

    h2p->Draw("colz");
    // TF1 *fpol = new TF1("fpol","[0]+[1]*x",0,4000);
    // h2p->Fit("fpol","WL");
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
