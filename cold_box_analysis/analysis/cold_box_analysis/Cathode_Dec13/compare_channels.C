#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void compare_channels(){

  vector<Int_t> channels = {1};
  TFile *f1 = new TFile("C1XArapuca_Efield_10kV_A4ch2_250MHz_cosmic00000/analyzed.root","READ");
  TFile *f2 = new TFile("C2XArapuca_Efield_10kV_A1ch1_250MHz_cosmic00000/analyzed.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  TTree *t2 = (TTree*)f2->Get("t1");
  vector<ADC_DATA> ch(channels.size());
  vector<ADC_DATA> ch2(channels.size());
  vector<TBranch*> bch(channels.size());
  vector<TBranch*> bch2(channels.size());

  vector<Double_t> a1;
  vector<Double_t> a2;
  for(Int_t k = 0; k<channels.size();k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
    bch2[k] = t2->GetBranch(Form("Ch%i",channels[k]));
    bch2[k]->SetAddress(&ch2[k]);
  }
  Int_t nentries = t1->GetEntries();
  for(Int_t i = 0; i<nentries; i++){
    for(Int_t k = 0; k<channels.size();k++){
      bch[k]->GetEvent(i);
      bch2[k]->GetEvent(i);
      a1.push_back(ch[k].charge);
      a2.push_back(ch2[k].charge);
  //     // code here
  //     // for(Int_t j = 0; j<memorydepth; j++){
  //     //   cout << ch[k].wvf[j] << endl;
  //     // }
    }
  }

  TGraph *g = new TGraph(a1.size(),&a1[0],&a2[0]);
  g->Draw("AP*");
}
