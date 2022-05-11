#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void averaged_waveform_spe(){
  string mychannel = "Ch2";
  vector<Double_t> volts = {5};
  // vector<Double_t> volts = {3,4,4.2,4.4,4.6,4.8,5,6,7.5};
  // vector<Double_t> volts = {10,12.5,15,17.5,20,22.5,25,27.5,30};
  // vector<Double_t> volts = {7.5};
  Int_t n = volts.size();

  vector<string> files(n);
  for(Int_t i=0; i<n; i++){
    if(static_cast<Int_t>(volts[i])-volts[i]!=0) files[i] = Form("%.1fV",volts[i]);
    else files[i] = Form("%.0fV",volts[i]);
  }
  
  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);

  TFile *fout = new TFile("few_pe_response.root","RECREATE");
  TH1D *hpe = new TH1D("hpe","hpe",memorydepth,0,4*memorydepth);

  vector<Double_t> avg(memorydepth);
  Double_t total_waveforms=0;
  for(Int_t i = 0; i<n; i++){
    files[i] = files[i]+conc;
    f[i] = new TFile(files[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch(mychannel.c_str());
    bch[i]->SetAddress(&ch[i]);


    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
    // for(Int_t j = 0; j<2000; j++){
      bch[i]->GetEvent(j);
      if(ch[i].selection==0 && ch[i].peak>60 && ch[i].peak<120){
        for(Int_t k = 0; k<memorydepth; k++){
          avg[k]+=ch[i].wvf[k];
        }
        total_waveforms+=1;
      }
    }
  }

  for(Int_t k = 0; k<memorydepth; k++){
    hpe->SetBinContent(k+1,avg[k]/total_waveforms);
  }



  fout->WriteTObject(hpe,"averaged_few_pe","TObject::kOverwrite");
  
  
}
