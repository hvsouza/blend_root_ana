#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void averaged_waveform_mu(){
  
  string mychannel = "Ch1";
  vector<Double_t> volts = {0};
  Int_t n = volts.size();

  vector<string> files(n);
  for(Int_t i=0; i<n; i++){
    if(static_cast<Int_t>(volts[i])-volts[i]!=0) files[i] = Form("%.1fV",volts[i]);
    else files[i] = Form("%.0fV",volts[i]);
  }
  
  Int_t nbins;
  Double_t dtime = 4;

  Int_t newStarting_point = 128;

  string conc = "analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);

  TFile *fout = new TFile("mu_response.root","RECREATE");
  TH1D *hmu = new TH1D("hmu","hmu",memorydepth-newStarting_point,0,4*memorydepth-4*newStarting_point);

  TH2D *hpersistance = new TH2D("hpersistance","hpersistance",memorydepth-newStarting_point,0,(memorydepth-newStarting_point)*4,2000,-2,2);
  vector<Double_t> avg(memorydepth-newStarting_point);
  Double_t total_waveforms=0;
  
  for(Int_t i = 0; i<n; i++){
    files[i] = conc;
    f[i] = new TFile(files[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch(mychannel.c_str());
    bch[i]->SetAddress(&ch[i]);


    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
    // for(Int_t j = 0; j<2000; j++){
      bch[i]->GetEvent(j);
      Bool_t triggered = false;
      // if(ch[i].selection==0 && ch[i].peak<1){
      if(ch[i].selection==0 && ch[i].peak < 0.8 && ch[i].peak>0.3){
        Int_t original_trigger = 0;
        for(Int_t k = 0; k<300; k++){
          if(ch[i].wvf[k]>=0.09 && triggered==false){
            triggered=true;
            original_trigger = k;
            // cout << j << " " << k << endl;
          }
        }
        Int_t shift = original_trigger - newStarting_point;
        for(Int_t k = 0; k<memorydepth-newStarting_point; k++){
          avg[k]+=ch[i].wvf[k+shift];
          hpersistance->Fill(k*4,ch[i].wvf[k+shift]);
        }
        total_waveforms+=1;
      }
    }
  }

  for(Int_t k = 0; k<memorydepth-newStarting_point; k++){
    // cout << k << " " << avg[k] << endl;
    hmu->SetBinContent(k+1,avg[k]/total_waveforms);
  }



  fout->WriteTObject(hmu,"averaged_mu","TObject::kOverwrite");
  fout->WriteTObject(hpersistance,"hpersistance","TObject::kOverwrite");
  
  
}
