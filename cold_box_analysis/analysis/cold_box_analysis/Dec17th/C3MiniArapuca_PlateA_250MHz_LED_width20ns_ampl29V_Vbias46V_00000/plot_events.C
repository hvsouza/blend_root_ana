#define memorydepth 938
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_events(){

  vector<Int_t> channels = {1};
  TFile *f1 = new TFile("analyzed.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  vector<ADC_DATA> ch(channels.size());
  vector<TBranch*> bch(channels.size());


  TH2D *h = new TH2D("h","h",938,0,938*4,10000,-0.4,0.4);
  ifstream fin;
  fin.open("valid_events.log");
  Int_t event = 0;
  for(Int_t k = 0; k<channels.size();k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
  }
  Int_t nentries = t1->GetEntries();
  while(!fin.fail()){
    fin >> event;
    if(fin.bad() || fin.fail()) break;
      
    for(Int_t k = 0; k<channels.size();k++){
      bch[k]->GetEvent(event);
      // code here
      for(Int_t j = 0; j<memorydepth; j++){
        h->Fill(j*4,ch[k].wvf[j]);
      }
    }
  }

  h->Draw("colz");
}
