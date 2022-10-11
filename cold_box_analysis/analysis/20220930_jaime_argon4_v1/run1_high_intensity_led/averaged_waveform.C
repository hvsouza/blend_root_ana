#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void averaged_waveform(){
  
  string mychannel = "Ch1";
  Int_t n = 1;

  vector<string> files(n);
  
  Int_t nbins;
  Double_t dtime = 4;



  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);

  TFile *fout = new TFile("response.root","RECREATE");
  TH1D *hmu = new TH1D("hmu","hmu",memorydepth,0,dtime*memorydepth);

  Double_t spe = 1400;
  vector<Double_t> avg(memorydepth);
  Double_t total_waveforms = 0;
  for(Int_t i = 0; i<n; i++){
    f[i] = new TFile("analyzed.root","READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch(mychannel.c_str());
    bch[i]->SetAddress(&ch[i]);

    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
    // for(Int_t j = 0; j<2000; j++){
      bch[i]->GetEvent(j);
      // if(ch[i].selection==0 && ch[i].peak<1){
      if(ch[i].selection==0 && ch[i].charge/spe>100 && ch[i].charge/spe < 120){
        for(Int_t k = 0; k<memorydepth; k++){
          avg[k]+=ch[i].wvf[k];
          }
        total_waveforms+=1;
        }
      }
    }

  for(Int_t k = 0; k<memorydepth; k++){
    // cout << k << " " << avg[k] << endl;
    hmu->SetBinContent(k+1,avg[k]/total_waveforms);
  }



  fout->WriteTObject(hmu,"averaged_mu","TObject::kOverwrite");

  }
