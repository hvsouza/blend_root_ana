const Int_t maxHits = 1328;

class DATA{
  Int_t subrun;
  Int_t event;
  Double_t evttime;
  Bool_t isdata;
  Double_t taulife;
  UInt_t triggernumber;
  Double_t triggertime;
  Int_t no_hits;
  Int_t no_hits_stored;
  Float_t hit_charge[maxHits];
  Float_t hit_ph[maxHits];
  Float_t hit_goodnessOfFit[maxHits];  




};

void read_tree(){



  TFile *f = new TFile("ana_hist.root","READ");
  TDirectoryFile *fd = (TDirectoryFile*)f->Get("analysistree");
  fd->cd();
  TTree *t1 = (TTree*)fd->Get("anatree");


}
