/* 
 To this values, do:
analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
 */
#define maxHits 1328
#define maxTracks 9


TFile *f;
TDirectoryFile *fd;
TTree *t1;


TBranch *b_subrun;
TBranch *b_event;
TBranch *b_evttime;
TBranch *b_isdata;
TBranch *b_taulife;
TBranch *b_triggernumber;
TBranch *b_triggertime;
TBranch *b_no_hits;
TBranch *b_no_hits_stored;
TBranch *b_hit_plane;
TBranch *b_hit_wire;
TBranch *b_hit_channel;
TBranch *b_hit_tpc;
TBranch *b_hit_charge;
TBranch *b_hit_ph;
TBranch *b_hit_goodnessOfFit;
TBranch *b_hit_trkid ;

TBranch *b_ntracks_pandoraTrack;
TBranch *b_trkId_pandoraTrack;
TBranch *b_trkdqdx_pandoraTrack;
TBranch *b_trkstartx_pandoraTrack;
TBranch *b_trkstarty_pandoraTrack;
TBranch *b_trkstartz_pandoraTrack;
TBranch *b_trkendx_pandoraTrack;
TBranch *b_trkendy_pandoraTrack;
TBranch *b_trkendz_pandoraTrack;
TBranch *b_trktheta_pandoraTrack;
TBranch *b_trkphi_pandoraTrack;


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
  Short_t hit_wire[maxHits];
  Short_t hit_channel[maxHits];
  Short_t hit_plane[maxHits];
  Short_t hit_tpc[maxHits];
  Float_t hit_charge[maxHits];
  Float_t hit_ph[maxHits];
  Float_t hit_goodnessOfFit[maxHits];  
  Short_t hit_trkid[maxHits];  

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

  void read_tree(){
    f = new TFile("ana_hist.root","READ");
    fd = (TDirectoryFile*)f->Get("analysistree");
    fd->cd();
    t1 = (TTree*)fd->Get("anatree");



    b_subrun = t1->GetBranch("subrun");
    b_event = t1->GetBranch("event");
    b_evttime = t1->GetBranch("evttime");
    b_isdata = t1->GetBranch("isdata");
    b_taulife = t1->GetBranch("taulife");
    b_triggernumber = t1->GetBranch("triggernumber");
    b_triggertime = t1->GetBranch("triggertime");
    b_no_hits = t1->GetBranch("no_hits");
    b_no_hits_stored = t1->GetBranch("no_hits_stored");
    b_hit_plane = t1->GetBranch("hit_plane");
    b_hit_wire = t1->GetBranch("hit_wire");
    b_hit_channel = t1->GetBranch("hit_channel");
    b_hit_tpc = t1->GetBranch("hit_tpc");
    b_hit_charge = t1->GetBranch("hit_charge");
    b_hit_ph = t1->GetBranch("hit_ph");
    b_hit_goodnessOfFit = t1->GetBranch("hit_goodnessOfFit");
    b_hit_trkid  = t1->GetBranch("hit_trkid");

    b_ntracks_pandoraTrack = t1->GetBranch("ntracks_pandoraTrack");
    b_trkId_pandoraTrack = t1->GetBranch("trkId_pandoraTrack");
    b_trkdqdx_pandoraTrack = t1->GetBranch("trkdqdx_pandoraTrack");
    b_trkstartx_pandoraTrack = t1->GetBranch("trkstartx_pandoraTrack");
    b_trkstarty_pandoraTrack = t1->GetBranch("trkstarty_pandoraTrack");
    b_trkstartz_pandoraTrack = t1->GetBranch("trkstartz_pandoraTrack");
    b_trkendx_pandoraTrack = t1->GetBranch("trkendx_pandoraTrack");
    b_trkendy_pandoraTrack = t1->GetBranch("trkendy_pandoraTrack");
    b_trkendz_pandoraTrack = t1->GetBranch("trkendz_pandoraTrack");
    b_trktheta_pandoraTrack = t1->GetBranch("trktheta_pandoraTrack");
    b_trkphi_pandoraTrack = t1->GetBranch("trkphi_pandoraTrack");

    b_subrun->SetAddress(&subrun);
    b_event->SetAddress(&event);
    b_evttime->SetAddress(&evttime);
    b_isdata->SetAddress(&isdata);
    b_taulife->SetAddress(&taulife);
    b_triggernumber->SetAddress(&triggernumber);
    b_triggertime->SetAddress(&triggertime);
    b_no_hits->SetAddress(&no_hits);
    b_no_hits_stored->SetAddress(&no_hits_stored);
    b_hit_wire->SetAddress(&hit_wire);
    b_hit_channel->SetAddress(&hit_channel);
    b_hit_plane->SetAddress(&hit_plane);
    b_hit_tpc->SetAddress(&hit_tpc);
    b_hit_charge->SetAddress(&hit_charge);
    b_hit_ph->SetAddress(&hit_ph);
    b_hit_goodnessOfFit->SetAddress(&hit_goodnessOfFit);
    b_hit_trkid->SetAddress(&hit_trkid);

    b_ntracks_pandoraTrack->SetAddress(&ntracks_pandoraTrack);
    b_trkId_pandoraTrack->SetAddress(&trkId_pandoraTrack);
    b_trkdqdx_pandoraTrack->SetAddress(&trkdqdx_pandoraTrack);
    b_trkstartx_pandoraTrack->SetAddress(&trkstartx_pandoraTrack);
    b_trkstarty_pandoraTrack->SetAddress(&trkstarty_pandoraTrack);
    b_trkstartz_pandoraTrack->SetAddress(&trkstartz_pandoraTrack);
    b_trkendx_pandoraTrack->SetAddress(&trkendx_pandoraTrack);
    b_trkendy_pandoraTrack->SetAddress(&trkendy_pandoraTrack);
    b_trkendz_pandoraTrack->SetAddress(&trkendz_pandoraTrack);
    b_trktheta_pandoraTrack->SetAddress(&trktheta_pandoraTrack);
    b_trkphi_pandoraTrack->SetAddress(&trkphi_pandoraTrack);
  }


  void plot_tracks(){
    // vector<Double_t> 
    
  }

  void test_some_events(Int_t mevent=2, Int_t mtpc = 2, Int_t mplane = 0){
    
    TH1D *h = new TH1D("h","h",200,0,0);
    vector<TH1D*> ht(3);
    for(Int_t i = 0; i<3; i++){
      ht[i] = new TH1D(Form("ht%d",i),Form("ht%d",i),200,0,0);
    }

    cout << "total_tracks current_track plane wire charge" << endl;
    for(Int_t l = 0; l<t1->GetEntries(); l++){
      t1->GetEntry(l);

      if(l==mevent){
      
        for(Int_t i = 0; i<no_hits_stored; i++){
          // printf("%.0f\n",evttime);
          if(hit_tpc[i]==mtpc){
            // This prints all the traks
            // for(Int_t g = 0; g<ntracks_pandoraTrack; g++){
            //   cout << trkId_pandoraTrack[g] << " ";
            // }
            cout << ntracks_pandoraTrack << " " << hit_trkid[i] << " " << hit_plane[i] << " " << hit_wire[i] << " " << hit_charge[i]  << endl;
            h->Fill(hit_charge[i]);
          }
        
        }
        // cout << "\n\n" << endl;
      }
      for(Int_t i = 0; i<ntracks_pandoraTrack; i++){
        // for(Int_t i = 0; i<1; i++){
        for(Int_t j = 0; j<3; j++){
          for(Int_t k = 0; k<2000; k++){
            // cout << trkId_pandoraTrack[i] << " " << j << " " << k << " " << trkdqdx_pandoraTrack[i][j][k] << endl;
            if(trkdqdx_pandoraTrack[i][j][k]!=0) ht[j]->Fill(trkdqdx_pandoraTrack[i][j][k]);
          }
        }
      }
    }
    // cout << event << endl;
    h->Draw();

    TCanvas *c2 = new TCanvas();
    THStack *hs = new THStack("hs","hs");
    hs->Add(ht[0]);
    hs->Add(ht[1]);
    hs->Add(ht[2]);

    hs->Draw("nostack");
  }


};




