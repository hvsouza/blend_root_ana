/* 
 To this values, do:
analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
 */
#define maxHits 1617
#define maxTracks 10


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
TBranch *b_hit_peakT;
TBranch *b_hit_tpc;
TBranch *b_hit_charge;
TBranch *b_hit_ph;
TBranch *b_hit_goodnessOfFit;
TBranch *b_hit_trkid ;

TBranch *b_ntracks_pandoraTrack;
TBranch *b_trkId_pandoraTrack;
TBranch *b_trkdqdx_pandoraTrack;
TBranch *b_ntrkhits_pandoraTrack;
TBranch *b_trkstartx_pandoraTrack;
TBranch *b_trkstarty_pandoraTrack;
TBranch *b_trkstartz_pandoraTrack;
TBranch *b_trkendx_pandoraTrack;
TBranch *b_trkendy_pandoraTrack;
TBranch *b_trkendz_pandoraTrack;
TBranch *b_trktheta_pandoraTrack;
TBranch *b_trkphi_pandoraTrack;

// Double_t xlenght = 337.5;
// Double_t ylenght = 337.5;
// Double_t zlenght = 337.5;

vector<Double_t> wallx = {-170,0}; //in cm // folling Laura Z. standard
vector<Double_t> wally = {0,300}; //in cm
vector<Double_t> wallz = {-10,25}; //in cm
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
  Float_t trkdqdx_pandoraTrack[maxTracks][3][2000];
  Float_t trkstartx_pandoraTrack[maxTracks];
  Float_t trkstarty_pandoraTrack[maxTracks];
  Float_t trkstartz_pandoraTrack[maxTracks];
  Float_t trkendx_pandoraTrack[maxTracks];
  Float_t trkendy_pandoraTrack[maxTracks];
  Float_t trkendz_pandoraTrack[maxTracks];
  Float_t trktheta_pandoraTrack[maxTracks];
  Float_t trkphi_pandoraTrack[maxTracks];


  Float_t *startz = trkstartx_pandoraTrack;
  Float_t *startx = trkstarty_pandoraTrack;
  Float_t *starty = trkstartz_pandoraTrack;
  Float_t *endz = trkendx_pandoraTrack;
  Float_t *endx = trkendy_pandoraTrack;
  Float_t *endy = trkendz_pandoraTrack;


  vector<Color_t> planesColor = {kRed,kBlue,kMagenta+2};
  vector<string> planesName{"U Induction","Y Induction","Z Collection"};
  
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
    b_hit_peakT = t1->GetBranch("hit_peakT");
    b_hit_tpc = t1->GetBranch("hit_tpc");
    b_hit_charge = t1->GetBranch("hit_charge");
    b_hit_ph = t1->GetBranch("hit_ph");
    b_hit_goodnessOfFit = t1->GetBranch("hit_goodnessOfFit");
    b_hit_trkid  = t1->GetBranch("hit_trkid");

    b_ntracks_pandoraTrack = t1->GetBranch("ntracks_pandoraTrack");
    b_trkId_pandoraTrack = t1->GetBranch("trkId_pandoraTrack");
    b_ntrkhits_pandoraTrack = t1->GetBranch("ntrkhits_pandoraTrack");
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
    b_hit_peakT->SetAddress(&hit_peakT);
    b_hit_plane->SetAddress(&hit_plane);
    b_hit_tpc->SetAddress(&hit_tpc);
    b_hit_charge->SetAddress(&hit_charge);
    b_hit_ph->SetAddress(&hit_ph);
    b_hit_goodnessOfFit->SetAddress(&hit_goodnessOfFit);
    b_hit_trkid->SetAddress(&hit_trkid);

    b_ntracks_pandoraTrack->SetAddress(&ntracks_pandoraTrack);
    b_trkId_pandoraTrack->SetAddress(&trkId_pandoraTrack);
    b_ntrkhits_pandoraTrack->SetAddress(&ntrkhits_pandoraTrack);
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


  void plot_tracks(Int_t mevent=-99){
    Bool_t not_selecting_event = false;
    Bool_t not_selecting_tpc = false;
    Width_t linewidth = 3;
    if(mevent==-99){
      not_selecting_event = true;
      linewidth = 1;
    }
    vector<Double_t> blankx = {wallx[0],wallx[1]};
    vector<Double_t> blanky = {wally[0],wally[1]};
    vector<Double_t> blankz = {wallz[0],wallz[1]};
    TCanvas *ccb = new TCanvas("ccb","ccb",1920,0,1920,680);
    TGraph2D *gphantom = new TGraph2D(2,&blankx[0],&blanky[0],&blankz[0]);
    gphantom->SetNameTitle("TDE","TDE");
    gphantom->GetXaxis()->SetTitle("x");
    gphantom->GetYaxis()->SetTitle("y");
    gphantom->GetZaxis()->SetTitle("z");
    gphantom->Draw("P");

    vector<TPolyLine3D*> tracks;
    Int_t ntracks=0;
    for(Int_t l = 0; l<t1->GetEntries(); l++){
      t1->GetEntry(l);

      if(event==mevent || not_selecting_event){

        for(Int_t itr = 0; itr<ntracks_pandoraTrack; itr++){
          Float_t tracklength = startz[itr] - endz[itr];//the values are negative, so watchout.
          endz[itr] = 23 - abs(tracklength); 
          startz[itr] = 23;// assuming always enter on the top

                  
          tracks.push_back(new TPolyLine3D(2));
          tracks[ntracks]->SetPoint(0,startx[itr], starty[itr], startz[itr]);
          tracks[ntracks]->SetPoint(1,endx[itr], endy[itr], endz[itr]);
          // cout << hit_tpc[i] << endl;
          tracks[ntracks]->SetLineWidth(linewidth);

          // if(selection_dqdx(itr)){
            tracks[ntracks]->SetLineColor(kBlue);      
          // }
          // else{
            // tracks[ntracks]->SetLineColor(kRed);      
          // }
          ntracks++;
          
        }

        
      }
      
    }

    for(Int_t k = 0; k<tracks.size(); k++){
      tracks[k]->Draw("P");
      // tracks[k]->Print();
    }


    return;


  }

void printInfo(string name, Float_t variable){
    cout << name << ": " << variable << " ";
  }
  void test_some_events(Int_t mevent=-99, Int_t mtpc = -99, Int_t mplane = -99){

    delete gROOT->FindObject("h");
    delete gROOT->FindObject("ht0");
    delete gROOT->FindObject("ht1");
    delete gROOT->FindObject("ht2");
    delete gROOT->FindObject("c1");
    delete gROOT->FindObject("c2");
    TH1D *h = new TH1D("h","h",200,0,0);
    TH1D *hdqdx = new TH1D("hdqdx","hdqdx",200,0,1600);
    vector<TH1D*> ht(3);
    for(Int_t i = 0; i<3; i++){
      ht[i] = new TH1D(Form("ht%d",i),Form("ht%d",i),200,0,1600);
      ht[i]->SetNameTitle(Form("dQ/dx %s",planesName[i].c_str()),Form("dQ/dx %s",planesName[i].c_str()));
      ht[i]->SetLineColor(planesColor[i]);
    }

    Bool_t not_selecting_event = false;
    Bool_t not_selecting_tpc = false;
    Bool_t not_selecting_plane = false;

    if(mevent==-99) not_selecting_event = true;
    if(mtpc==-99) not_selecting_tpc = true;
    if(mplane==-99) not_selecting_plane = true;


    // cout << "hit_trkid hit_plane hit_wire hit_peakT" << endl;
    Int_t starting_event = 0;
    Int_t finish_event = t1->GetEntries();
    if(not_selecting_event==false)
      {
        starting_event = mevent; 
        finish_event = mevent+1; 
      }
    for(Int_t l =starting_event; l<finish_event; l++){
      t1->GetEntry(l);


      Int_t trackRef = -99;

      vector<Bool_t> passTPC1={false,false,false};
      vector<Bool_t> passTPC2={false,false,false};
    
      for(Int_t i = 0; i<no_hits_stored; i++){
        // printf("%.2f\n",evttime);
        // checkBothTPCsHit(passTPC1,passTPC2,i);
        
        if(hit_tpc[i]==mtpc || not_selecting_tpc){
          // This prints all the tracks
          // for(Int_t g = 0; g<ntracks_pandoraTrack; g++){
          //   cout << trkId_pandoraTrack[g] << " ";
          // }
          // if(hit_plane[i]==0) cout << hit_trkid[i] << " " << hit_plane[i] << " " << hit_wire[i] << " " << hit_peakT[i]  << endl;
          if(hit_plane[i]==mplane || not_selecting_plane)
            {
              // cout << trackRef << " " << hit_trkid[i] << endl;
              if(trackRef!=hit_trkid[i]){               
                // trackRef = hit_trkid[i];
                // printInfo("event",event);
                // printInfo("Id",hit_trkid[i]);
                // printInfo("wire",hit_wire[i]);
                // printInfo("peakT",hit_peakT[i]);
                // printInfo("charge",hit_charge[i]);
                // cout << "\n";
              }

            }
          h->Fill(hit_charge[i]);

        }
        
      }



      for(Int_t i2 = 0; i2<ntracks_pandoraTrack; i2++){
        // for(Int_t i = 0; i<1; i++){

        Double_t trackTotaldqdx = 0;
        for(Int_t j2 = 0; j2<3; j2++){          
          for(Int_t k2 = 0; k2<ntrkhits_pandoraTrack[i2][j2]; k2++){

            ht[j2]->Fill(trkdqdx_pandoraTrack[i2][j2][k2]);
            cout << trkId_pandoraTrack[i2] << " " << j2 << " " << k2 << " " << trkdqdx_pandoraTrack[i2][j2][k2] << " " << ntrkhits_pandoraTrack[i2][j2] << endl;
              // trackTotaldqdx += trkdqdx_pandoraTrack[i2][j2][k2];
          }
        }

      }

      // break;
          
      // cout << "\n" << endl;
      
      
    }
    // cout << "\n" << endl;
    TCanvas *c1 = new TCanvas("c1","c1");
    h->Draw();

    TCanvas *c2 = new TCanvas("c2","c2");
    THStack *hs = new THStack("hs","hs");
    
    hs->Add(ht[0]);
    hs->Add(ht[1]);
    hs->Add(ht[2]);

    hs->Draw("nostack");
    c2->BuildLegend();
    // delete h;
    // delete ht[0];
    // delete ht[1];
    // delete ht[2];

    TCanvas *ctdqdx = new TCanvas("ctdqdx","ctdqdx");
    hdqdx->Add(ht[0]);
    hdqdx->Add(ht[1]);
    hdqdx->Add(ht[2]);
    hdqdx->Draw();

  }


  // Bool_t selection_dqdx(Int_t i2){
  //   vector<Double_t> planes_threshold= {200,300,300};
  //   for(Int_t j2 = 0; j2<3; j2++){
  //     for(Int_t k2 = 0; k2<2000; k2++){
  //       if(trkdqdx_pandoraTrack[i2][j2][k2]!=0){
  //         // trkdqdx_pandoraTrack[i2][j2][k2];
  //         // trackTotaldqdx += trkdqdx_pandoraTrack[i2][j2][k2];
  //       }
  //       else break;
  //     }
  //   }
   
  // }
  
  void checkBothTPCsHit(vector<Bool_t> &passTPC1, vector<Bool_t> &passTPC2, Int_t i)
  {
      
    if(hit_trkid[i]==0){
      if(hit_tpc[i]==0){
        for(Int_t idp = 0; idp<3; idp++){
          if(hit_plane[i]==idp) passTPC1[idp] = true;
        }
      }
      if(hit_tpc[i]==2){
        for(Int_t idp = 0; idp<3; idp++){
          if(hit_plane[i]==idp) passTPC2[idp] = true;
        }
      }
      for(Int_t idp = 0; idp<3; idp++){
        if(passTPC1[idp] == true && passTPC2[idp] == true){
          cout << "WARNING!!!! \n\n\n " << endl;
          printInfo("event",event);
          printInfo("plane",idp);
          printInfo("Id",hit_trkid[i]);
          printInfo("wire",hit_wire[i]);
          printInfo("peakT",hit_peakT[i]);
          cout << "\n";
        }
        // else{
          // cout << passTPC1[idp] << " " << passTPC2[idp] << " ";
        // }
      }
      // cout << "\n";
    }
  }


};




