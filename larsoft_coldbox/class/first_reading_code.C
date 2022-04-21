/* 
   To this values, do:
   analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
*/

const Int_t maxHits = 1617;
const Int_t maxTracks = 10;

class DATA{

public:
  Int_t subrun;
  Int_t event;
  Double_t evttime; //epoch time in seconds
  Char_t isdata;
  Double_t taulife;
  UInt_t triggernumber;
  Double_t triggertime;

  Int_t no_hits;
  Int_t no_hits_stored;
  Short_t hit_wire[maxHits];
  Short_t hit_channel[maxHits];
  Float_t hit_peakT[maxHits];
  Short_t hit_plane[maxHits];
  Short_t hit_tpc[maxHits];
  Float_t hit_charge[maxHits];
  Float_t hit_ph[maxHits];
  Float_t hit_goodnessOfFit[maxHits];  
  Short_t hit_trkid[maxHits];  

  Short_t ntracks_pandoraTrack;
  Short_t trkId_pandoraTrack[maxTracks];
  Short_t ntrkhits_pandoraTrack[maxTracks][3];
  Short_t trkmomrange_pandoraTrack[maxTracks];
  Float_t trkdqdx_pandoraTrack[maxTracks][3][2000];
  Float_t trkxyz_pandoraTrack[maxTracks][3][2000][3];
  Float_t trkstartx_pandoraTrack[maxTracks];
  Float_t trkstarty_pandoraTrack[maxTracks];
  Float_t trkstartz_pandoraTrack[maxTracks];
  Float_t trkendx_pandoraTrack[maxTracks];
  Float_t trkendy_pandoraTrack[maxTracks];
  Float_t trkendz_pandoraTrack[maxTracks];
  Float_t trktheta_pandoraTrack[maxTracks];
  Float_t trkphi_pandoraTrack[maxTracks];

  TFile *f;
  TDirectoryFile *fd;
  TTree *t1;


  void read_tree(){
    f = new TFile("ana_hist_429_1.root","READ");
    fd = (TDirectoryFile*)f->Get("analysistree");
    fd->cd();
    t1 = (TTree*)fd->Get("anatree");


    t1->SetBranchAddress("subrun",&subrun);
    t1->SetBranchAddress("event",&event);
    t1->SetBranchAddress("evttime",&evttime);
    t1->SetBranchAddress("isdata",&isdata);
    t1->SetBranchAddress("taulife",&taulife);
    t1->SetBranchAddress("triggernumber",&triggernumber);
    t1->SetBranchAddress("triggertime",&triggertime);
    t1->SetBranchAddress("no_hits",&no_hits);
    t1->SetBranchAddress("no_hits_stored",&no_hits_stored);
    t1->SetBranchAddress("hit_plane",&hit_plane);
    t1->SetBranchAddress("hit_wire",&hit_wire);
    t1->SetBranchAddress("hit_channel",&hit_channel);
    t1->SetBranchAddress("hit_peakT",&hit_peakT);
    t1->SetBranchAddress("hit_tpc",&hit_tpc);
    t1->SetBranchAddress("hit_charge",&hit_charge);
    t1->SetBranchAddress("hit_ph",&hit_ph);
    t1->SetBranchAddress("hit_goodnessOfFit",&hit_goodnessOfFit);
    t1->SetBranchAddress("hit_trkid",&hit_trkid);

    t1->SetBranchAddress("ntracks_pandoraTrack",&ntracks_pandoraTrack);
    t1->SetBranchAddress("trkId_pandoraTrack",&trkId_pandoraTrack);
    t1->SetBranchAddress("ntrkhits_pandoraTrack",&ntrkhits_pandoraTrack);
    t1->SetBranchAddress("trkmomrange_pandoraTrack",&trkmomrange_pandoraTrack);
    t1->SetBranchAddress("trkdqdx_pandoraTrack",&trkdqdx_pandoraTrack);
    t1->SetBranchAddress("trkxyz_pandoraTrack",&trkxyz_pandoraTrack);
    t1->SetBranchAddress("trkstartx_pandoraTrack",&trkstartx_pandoraTrack);
    t1->SetBranchAddress("trkstarty_pandoraTrack",&trkstarty_pandoraTrack);
    t1->SetBranchAddress("trkstartz_pandoraTrack",&trkstartz_pandoraTrack);
    t1->SetBranchAddress("trkendx_pandoraTrack",&trkendx_pandoraTrack);
    t1->SetBranchAddress("trkendy_pandoraTrack",&trkendy_pandoraTrack);
    t1->SetBranchAddress("trkendz_pandoraTrack",&trkendz_pandoraTrack);
    t1->SetBranchAddress("trktheta_pandoraTrack",&trktheta_pandoraTrack);
    t1->SetBranchAddress("trkphi_pandoraTrack",&trkphi_pandoraTrack);



    TH1D *h = new TH1D("h","h",200,0,0);
    vector<TH1D*> ht(3);
  
    for(Int_t i = 0; i<3; i++){
      ht[i] = new TH1D(Form("ht%d",i),Form("ht%d",i),200,0,0);
    }
  
    for(Int_t l = 0; l<t1->GetEntries(); l++){
      t1->GetEntry(l);    
      for(Int_t i = 0; i<no_hits_stored; i++){

        if(hit_tpc[i]==2){
          for(Int_t g = 0; g<ntracks_pandoraTrack; g++){
            cout << trkId_pandoraTrack[g] << " ";
          }
          cout << ntracks_pandoraTrack << " " << hit_trkid[i] << " " << hit_plane[i] << " " << hit_wire[i] <<  " " << hit_channel[i] << " " << hit_charge[i]  << endl;
          h->Fill(hit_charge[i]);
        }
        
      }
      cout << "\n\n" << endl;
    
      for(Int_t i = 0; i<ntracks_pandoraTrack; i++){
        for(Int_t j = 0; j<3; j++){
          for(Int_t k = 0; k<ntrkhits_pandoraTrack[i][j]; k++){
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

void first_reading_code(){
  DATA data;

  data.read_tree();
}


