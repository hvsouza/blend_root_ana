const Int_t maxHits = 1328;
const Int_t maxTracks = 9;

class DATA{

public:
  Int_t subrun;
  Int_t event;
  Double_t evttime; //epoch time in seconds
  Bool_t isdata;
  Double_t taulife;
  UInt_t triggernumber;
  Double_t triggertime;
  Int_t no_hits;
  Int_t no_hits_stored;
  Float_t hit_charge[maxHits];
  Float_t hit_ph[maxHits];
  Float_t hit_goodnessOfFit[maxHits];  

  Short_t ntracks_pandoraTrack;
  Short_t trkId_pandoraTrack[maxTracks];
  Float_t trkdqdx_pandoraTrack[maxTracks][3][2000];
  Float_t trkstartx_pandoraTrack[maxTracks];
  Float_t trkstarty_pandoraTrack[maxTracks];
  Float_t trkstartz_pandoraTrack[maxTracks];
  Float_t trkendx_pandoraTrack[maxTracks];
  Float_t trkendy_pandoraTrack[maxTracks];
  Float_t trkendz_pandoraTrack[maxTracks];
  Float_t trktheta_pandoraTrack[maxTracks];
  Float_t trkphi_pandoraTrack[maxTracks];



};

void read_tree(){

  TFile *f = new TFile("ana_hist.root","READ");
  TDirectoryFile *fd = (TDirectoryFile*)f->Get("analysistree");
  fd->cd();
  TTree *t1 = (TTree*)fd->Get("anatree");

  DATA data;

  TBranch *b_subrun = t1->GetBranch("subrun");
  TBranch *b_event = t1->GetBranch("event");
  TBranch *b_evttime = t1->GetBranch("evttime");
  TBranch *b_isdata = t1->GetBranch("isdata");
  TBranch *b_taulife = t1->GetBranch("taulife");
  TBranch *b_triggernumber = t1->GetBranch("triggernumber");
  TBranch *b_triggertime = t1->GetBranch("triggertime");
  TBranch *b_no_hits = t1->GetBranch("no_hits");
  TBranch *b_no_hits_stored = t1->GetBranch("no_hits_stored");
  TBranch *b_hit_charge = t1->GetBranch("hit_charge");
  TBranch *b_hit_ph = t1->GetBranch("hit_ph");
  TBranch *b_hit_goodnessOfFit  = t1->GetBranch("hit_goodnessOfFit");

  TBranch *b_ntracks_pandoraTrack = t1->GetBranch("ntracks_pandoraTrack");
  TBranch *b_trkId_pandoraTrack = t1->GetBranch("trkId_pandoraTrack");
  TBranch *b_trkdqdx_pandoraTrack = t1->GetBranch("trkdqdx_pandoraTrack");
  TBranch *b_trkstartx_pandoraTrack = t1->GetBranch("trkstartx_pandoraTrack");
  TBranch *b_trkstarty_pandoraTrack = t1->GetBranch("trkstarty_pandoraTrack");
  TBranch *b_trkstartz_pandoraTrack = t1->GetBranch("trkstartz_pandoraTrack");
  TBranch *b_trkendx_pandoraTrack = t1->GetBranch("trkendx_pandoraTrack");
  TBranch *b_trkendy_pandoraTrack = t1->GetBranch("trkendy_pandoraTrack");
  TBranch *b_trkendz_pandoraTrack = t1->GetBranch("trkendz_pandoraTrack");
  TBranch *b_trktheta_pandoraTrack = t1->GetBranch("trktheta_pandoraTrack");
  TBranch *b_trkphi_pandoraTrack = t1->GetBranch("trkphi_pandoraTrack");

  b_subrun->SetAddress(&data.subrun);
  b_event->SetAddress(&data.event);
  b_evttime->SetAddress(&data.evttime);
  b_isdata->SetAddress(&data.isdata);
  b_taulife->SetAddress(&data.taulife);
  b_triggernumber->SetAddress(&data.triggernumber);
  b_triggertime->SetAddress(&data.triggertime);
  b_no_hits->SetAddress(&data.no_hits);
  b_no_hits_stored->SetAddress(&data.no_hits_stored);
  b_hit_charge->SetAddress(&data.hit_charge);
  b_hit_ph->SetAddress(&data.hit_ph);
  b_hit_goodnessOfFit->SetAddress(&data.hit_goodnessOfFit);

  b_ntracks_pandoraTrack->SetAddress(&data.ntracks_pandoraTrack);
  b_trkId_pandoraTrack->SetAddress(&data.trkId_pandoraTrack);
  b_trkdqdx_pandoraTrack->SetAddress(&data.trkdqdx_pandoraTrack);
  b_trkstartx_pandoraTrack->SetAddress(&data.trkstartx_pandoraTrack);
  b_trkstarty_pandoraTrack->SetAddress(&data.trkstarty_pandoraTrack);
  b_trkstartz_pandoraTrack->SetAddress(&data.trkstartz_pandoraTrack);
  b_trkendx_pandoraTrack->SetAddress(&data.trkendx_pandoraTrack);
  b_trkendy_pandoraTrack->SetAddress(&data.trkendy_pandoraTrack);
  b_trkendz_pandoraTrack->SetAddress(&data.trkendz_pandoraTrack);
  b_trktheta_pandoraTrack->SetAddress(&data.trktheta_pandoraTrack);
  b_trkphi_pandoraTrack->SetAddress(&data.trkphi_pandoraTrack);


  TH1D *h = new TH1D("h","h",200,0,0);
  vector<TH1D*> ht(3);
  for(Int_t i = 0; i<3; i++){
    ht[i] = new TH1D(Form("ht%d",i),Form("ht%d",i),200,0,0);
  }
  
  for(Int_t l = 0; l<t1->GetEntries(); l++){
    t1->GetEntry(l);
    
    for(Int_t i = 0; i<data.no_hits_stored; i++){
      // printf("%.0f\n",data.evttime);
      h->Fill(data.hit_charge[i]);
    }
    for(Int_t i = 0; i<data.ntracks_pandoraTrack; i++){
    // for(Int_t i = 0; i<1; i++){
      for(Int_t j = 0; j<3; j++){
        for(Int_t k = 0; k<2000; k++){
          // cout << data.trkId_pandoraTrack[i] << " " << j << " " << k << " " << data.trkdqdx_pandoraTrack[i][j][k] << endl;
          if(data.trkdqdx_pandoraTrack[i][j][k]!=0) ht[j]->Fill(data.trkdqdx_pandoraTrack[i][j][k]);
        }
      }
    }
  }
  cout << data.event << endl;
  h->Draw();

  TCanvas *c2 = new TCanvas();
  THStack *hs = new THStack("hs","hs");
  hs->Add(ht[0]);
  hs->Add(ht[1]);
  hs->Add(ht[2]);

  hs->Draw("nostack");
}


