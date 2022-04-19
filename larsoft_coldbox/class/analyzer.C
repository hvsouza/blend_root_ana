/* 
   To this values, do:
   analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
*/


#include <fstream>
#include <iostream>
#include "Riostream.h"
#include <time.h>       // time_t, struct tm, difftime, time, mktime
#include <cmath>
#include <numeric>



#include "TROOT.h"
#include "TLatex.h"
#include <TMinuit.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TLine.h>
#include <TTimeStamp.h>


class DATA{

public:

  Int_t n_events;
  
  string file_name = "ana_hist_something.root";
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


  Float_t startz[maxTracks];
  Float_t startx[maxTracks];
  Float_t starty[maxTracks];
  Float_t endz[maxTracks];
  Float_t endx[maxTracks];
  Float_t endy[maxTracks];
  Float_t total_range[maxTracks];

  vector<vector<Bool_t>> selected_dqdx;
  vector<vector<Bool_t>> selected_phi;


  vector<Color_t> planesColor = {kRed,kBlue,kMagenta+2};
  vector<string> planesName{"U Induction","Y Induction","Z Collection"};

  vector<Double_t> wallx = {-170,0}; //in cm // folling Laura Z. standard
  vector<Double_t> wally = {0,300}; //in cm
  vector<Double_t> wallz = {-5,25}; //in cm
  Double_t pi = TMath::Pi();

  vector<vector<Double_t>> phi;
  vector<vector<Double_t>> theta;
  
  TFile *f;
  TDirectoryFile *fd;
  TTree *t1;

  Double_t phi_rotation_north_degree = 41.6; // north is at 41.6 degrees more or less 
  Double_t phi_rotation_north_rad = 41.6*pi/180;  


  DATA(string fname, Int_t nlimit);
  void read_tree();
  void FillHistoWithData(DATA &data);
  Double_t convert_north(Double_t mphi);


  Int_t nbins_htheta = 9;
  TH1D *hphi = new TH1D("hphi","hphi",36,0,2*pi);
  TH1D *hphi_north = new TH1D("hphi_north","hphi_north",32,0,2*pi);
  TH1D *htheta_t = new TH1D("htheta_t","htheta_t",nbins_htheta,0,90);
  TH1D *htheta = new TH1D("htheta","htheta",nbins_htheta,0,90);
  TH2D *hrange_theta = new TH2D("hrange_theta","hrange_theta",180,90,180,330,-10,320);

  Int_t maxPlotEvents = 500;
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
    TCanvas *ccb = new TCanvas("ccb","ccb",1920,0,1920,640);
    TGraph2D *gphantom = new TGraph2D(2,&blankx[0],&blanky[0],&blankz[0]);
    gphantom->SetNameTitle("TDE","TDE");
    gphantom->GetXaxis()->SetTitle("y (cm)");
    gphantom->GetYaxis()->SetTitle("z (cm)");
    gphantom->GetZaxis()->SetTitle("x (cm)");
    gphantom->Draw("P");

    Int_t nbinsx = abs(blankx[1]-blankx[0])*10;
    Int_t nbinsy = abs(blanky[1]-blanky[0])*10;
    Int_t nbinsz = abs(blankz[1]-blankz[0])*10;
    
    TH2D *hxy = new TH2D("hxy","hxy",nbinsx,blankx[0],blankx[1],nbinsy,blanky[0],blanky[1]);
    TH2D *hxz = new TH2D("hxz","hxz",nbinsx,blankx[0],blankx[1],nbinsz,blankz[0],blankz[1]);
    


    
    vector<TPolyLine3D*> tracks;
    Int_t ntracks=0;
    if(maxPlotEvents==0) maxPlotEvents = t1->GetEntries(); 
    for(Int_t l = 0; l<maxPlotEvents; l++){
      t1->GetEntry(l);

      if(event==mevent || not_selecting_event){

        for(Int_t itr = 0; itr<ntracks_pandoraTrack; itr++){
          getCoordinates(itr);
          Float_t tracklength = startz[itr] - endz[itr];//the values are negative, so watchout.
          // if(abs(tracklength)>26.5) cout << event << " " << trkId_pandoraTrack[itr] << " " <<  tracklength << " " << startz[itr] << " " << endz[itr] <<  endl;
          // cout << event << " " << trkId_pandoraTrack[itr] << " " << total_range[itr] << endl;
          
          if(selected_phi[l][itr]){
            for(Int_t j2 = 0; j2<3; j2++){
              for(Int_t k2 = 0; k2<ntrkhits_pandoraTrack[itr][j2]; k2++){
                Float_t zpoint = 23 - abs(trkxyz_pandoraTrack[itr][j2][k2][0]-startz[itr]);
                // cout << zpoint << endl;
                hxy->Fill(trkxyz_pandoraTrack[itr][j2][k2][1],trkxyz_pandoraTrack[itr][j2][k2][2]);
                hxz->Fill(trkxyz_pandoraTrack[itr][j2][k2][1],zpoint);
              }
            }
          }
  

          endz[itr] = 23 - abs(tracklength); 
          startz[itr] = 23;// assuming always enter on the top

          tracks.push_back(new TPolyLine3D(2));
          tracks[ntracks]->SetPoint(0,startx[itr], starty[itr], startz[itr]);
          tracks[ntracks]->SetPoint(1,endx[itr], endy[itr], endz[itr]);
          // cout << hit_tpc[i] << endl;
          tracks[ntracks]->SetLineWidth(linewidth);
          
          if(selected_phi[l][itr]){
            tracks[ntracks]->SetLineColor(kBlue);
          }
          else{
            tracks[ntracks]->SetLineColor(kRed);      
          }
          ntracks++;
          
          
        }

        
      }
      
    }

    for(Int_t k = 0; k<tracks.size(); k++){
      tracks[k]->Draw("P");
      // tracks[k]->Print();
    }

    TCanvas *cxy = new TCanvas("cxy","cxy",1920,650,960,386);
    cxy->cd();
    hxy->GetXaxis()->SetTitle("y (cm)");
    hxy->GetYaxis()->SetTitle("z (cm)");
    hxy->Draw("colz");

    TCanvas *cxz = new TCanvas("cxz","cxz",2880,650,960,386);
    cxz->cd();
    hxz->GetXaxis()->SetTitle("y (cm)");
    hxz->GetYaxis()->SetTitle("x (cm)");
    hxz->Draw("colz");
       

    TCanvas *crange_theta = new TCanvas("crange_theta","crange_theta");
    hrange_theta->Draw("colz");

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
    Int_t finish_event = n_events;
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
        vector<Double_t> total_dqdx = {0,0,0};
        for(Int_t j2 = 0; j2<3; j2++){          
          for(Int_t k2 = 0; k2<ntrkhits_pandoraTrack[i2][j2]; k2++){
            total_dqdx[j2]+=trkdqdx_pandoraTrack[i2][j2][k2];
            ht[j2]->Fill(trkdqdx_pandoraTrack[i2][j2][k2]);
            // cout << trkId_pandoraTrack[i2] << " " << j2 << " " << k2 << " " << trkdqdx_pandoraTrack[i2][j2][k2] << " " << ntrkhits_pandoraTrack[i2][j2] << endl;
            // trackTotaldqdx += trkdqdx_pandoraTrack[i2][j2][k2];
          }
          // ht[j2]->Fill(total_dqdx[j2]/ntrkhits_pandoraTrack[i2][j2]);
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
    
    // TCanvas *ctdqdx = new TCanvas("ctdqdx","ctdqdx");
    // hdqdx->Add(ht[0]);
    // hdqdx->Add(ht[1]);
    // hdqdx->Add(ht[2]);
    // hdqdx->Draw();

    TCanvas *cphi = new TCanvas("cphi","cphi");
    TAxis* aphi= hphi->GetXaxis();
    aphi->SetNdivisions(-504);
    aphi->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
    aphi->ChangeLabel(2,-1,-1,-1,-1,-1,"#pi/2");
    aphi->ChangeLabel(3,-1,-1,-1,-1,-1,"#pi");
    aphi->ChangeLabel(4,-1,-1,-1,-1,-1,"3#pi/2");
    aphi->ChangeLabel(-1,-1,-1,-1,-1,-1,"2#pi");
       
    hphi->Draw();

    TCanvas *cphi_north = new TCanvas("cphi_north","cphi_north");
    TAxis* aphi_north= hphi_north->GetXaxis();
    aphi_north->SetNdivisions(-504);
    aphi_north->ChangeLabel(1,-1,-1,-1,-1,-1,"N");
    aphi_north->ChangeLabel(2,-1,-1,-1,-1,-1,"W");
    aphi_north->ChangeLabel(3,-1,-1,-1,-1,-1,"S");
    aphi_north->ChangeLabel(4,-1,-1,-1,-1,-1,"E");
    aphi_north->ChangeLabel(-1,-1,-1,-1,-1,-1,"N");
       
    hphi_north->Draw();

    TCanvas *ctheta = new TCanvas("ctheta","ctheta");
    htheta_t->Scale(1./htheta->Integral("width"));
    htheta->Scale(1./htheta->Integral("width"));
    htheta_t->SetLineColor(kBlack);
    htheta->SetLineColor(kBlue);
    htheta_t->SetLineWidth(2);
    htheta->SetLineWidth(2);
    TAxis* a = htheta->GetXaxis();
    // a->SetNdivisions(-502);
    // a->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
    // a->ChangeLabel(2,-1,-1,-1,-1,-1,"#pi/4");
    // a->ChangeLabel(-1,-1,-1,-1,-1,-1,"#pi/2");
    htheta->GetYaxis()->SetTitle("dN/d#theta (degree^{-1})");
    htheta->GetXaxis()->SetTitle("#theta (degree)");
    htheta->Draw("");
    htheta_t->Draw("HIST SAME");

    // delete h;
    // delete ht[0];
    // delete ht[1];
    // delete ht[2];

  }


  void selection_dqdx(){

    for(Int_t l = 0; l<n_events; l++){
      t1->GetEntry(l);
      vector<Double_t> planes_threshold= {120,120,120};
      vector<Bool_t> selected = {false,false,false};
      Int_t ntrues=0;
      for(Int_t i2 = 0; i2<ntracks_pandoraTrack; i2++){
        for(Int_t j2 = 0; j2<3; j2++){
          for(Int_t k2 = 0; k2<ntrkhits_pandoraTrack[i2][j2]; k2++){
            if(trkdqdx_pandoraTrack[i2][j2][k2]>planes_threshold[j2]){
              if(selected[j2]==false){
                selected[j2] = true;
                ntrues++;
              }
            }
          }
        }
        if(ntrues>=2){
          selected_dqdx[l][i2] = true;
        }
        else{
          selected_dqdx[l][i2] = false;
        }
      }

    }
  }


  void getNewPhiTheta(Int_t l,Int_t i2){
    Double_t old_phi = trkphi_pandoraTrack[i2];
    if(old_phi<0) old_phi = 2*pi+old_phi;
    Double_t old_theta = trktheta_pandoraTrack[i2];

    
    Double_t xp = sin(old_phi)*sin(old_theta); //x' = y
    Double_t yp = cos(old_theta); // y' = z
    Double_t zp = cos(old_phi)*sin(old_theta); // z' = x
    Double_t nphi;
    Double_t ntheta;
    
    if(zp!=0) ntheta = atan(sqrt(xp*xp+yp*yp)/zp);
    else ntheta = pi/2;
    if(xp!=0)nphi = atan(yp/xp);
    else nphi = pi/2;
    if(ntheta<0){
      ntheta = pi+ntheta; // remember theta is negative here
      // nphi = nphi+pi;
    }
    
    // if(xp>=0 && yp>=0){ // this is not needed
      // cout << "good.... " << endl;
    // }
    // just to keep separate quadrants:
    if(xp<0 && yp>=0){ // second quadrant which gives negative phi
      nphi = nphi+pi;
    }
    else if(xp<0 && yp<=0){ // third quadrant which gives positive phi
      nphi = nphi+pi;
    }
    else if(xp>0 && yp<0){
      nphi = 2*pi+nphi;
    }
    if(ntheta<pi/2){
      ntheta = pi-ntheta;
      if(nphi<=pi) nphi = nphi + pi;
      else if(nphi>pi) nphi = nphi - pi;
      else{
        cout << "Something is very wrong... " << endl;
        return;
      }

      
    }

 
    // hphi->Fill(in_angle(nphi));
    // htheta->Fill(in_angle(ntheta));
    phi[l][i2] = nphi;
    theta[l][i2] = ntheta;
    
  }
  
  void selection_phi(Int_t cut = 1){
    
    for(Int_t l = 0; l<n_events; l++){
      t1->GetEntry(l);
      for(Int_t i2 = 0; i2<ntracks_pandoraTrack; i2++){
        selected_phi[l][i2] = false;
        Double_t mphi = in_angle(phi[l][i2]);
        switch(cut){
        case 1:
          if(mphi<=45 || mphi>=315) selected_phi[l][i2] = true;
          break;
        case 2:
          if(mphi<=90 || mphi>=270) selected_phi[l][i2] = true;
          break;
        }
      }   
    }
  }
  

  void change_pos(Float_t  &a,Float_t &b){
    Float_t tmp = a;
    a = b;
    b = tmp;
  }
  void getCoordinates(Int_t itr){

    total_range[itr] = 0;
    startz[itr] = trkstartx_pandoraTrack[itr];
    startx[itr] = trkstarty_pandoraTrack[itr];
    starty[itr] = trkstartz_pandoraTrack[itr];
    endz[itr] = trkendx_pandoraTrack[itr];
    endx[itr] = trkendy_pandoraTrack[itr];
    endy[itr] = trkendz_pandoraTrack[itr];
    if(startz[itr]>endz[itr]){
      change_pos(startz[itr],endz[itr]);
      change_pos(starty[itr],endy[itr]);
      change_pos(startx[itr],endx[itr]);
    }
    total_range[itr] += (startz[itr]-endz[itr])*(startz[itr]-endz[itr]);
    total_range[itr] += (starty[itr]-endy[itr])*(starty[itr]-endy[itr]);
    total_range[itr] += (startx[itr]-endx[itr])*(startx[itr]-endx[itr]);
    total_range[itr] = sqrt(total_range[itr]);
    
  }

  Double_t in_angle(Double_t val){
    return val*180/pi;
  }
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


DATA::DATA(string fname, Int_t nlimit = 0){
  
  file_name = fname;
  read_tree();
  if(nlimit==0) n_events = t1->GetEntries();
  else{
    n_events = nlimit;
    n_events = maxPlotEvents;
  }
  selected_phi.resize(n_events);
  selected_dqdx.resize(n_events);
  //getting phi and theta here is faster, I can do it only once
  phi.resize(n_events);
  theta.resize(n_events);
  
  for(Int_t l = 0; l<n_events; l++){
    t1->GetEntry(l);
    selected_phi[l].resize(ntracks_pandoraTrack,true);
    selected_dqdx[l].resize(ntracks_pandoraTrack,true);
    phi[l].resize(ntracks_pandoraTrack,0.);
    theta[l].resize(ntracks_pandoraTrack,0.);
    for(Int_t i2 = 0; i2<ntracks_pandoraTrack; i2++) getNewPhiTheta(l,i2);
  }
  
}


void DATA::FillHistoWithData(DATA &data){

  Double_t safe_distance = 5; // 20 cm away from the walls
  // data.selection_phi();
  for(Int_t l = 0; l<data.n_events; l++){
    data.t1->GetEntry(l);
    for(Int_t i2 = 0; i2<data.ntracks_pandoraTrack; i2++){
      data.getCoordinates(i2);
      Double_t mtheta = data.theta[l][i2];
      Double_t mphi = data.phi[l][i2];
      hrange_theta->Fill(in_angle(mtheta),data.total_range[i2]);
      htheta_t->Fill(in_angle(pi-(mtheta)));
      Bool_t hist_filled = false;
      if(data.total_range[i2]>15){
        // if(data.startx[i2]>(wallx[0]+safe_distance) && data.startx[i2]<(wallx[1]-safe_distance)){
        // if(data.starty[i2]>(wally[0]+safe_distance) && data.starty[i2]<(wally[1]-safe_distance)){
        // if(in_angle(mtheta)<120){
          hphi->Fill(mphi);
          hphi_north->Fill(convert_north(mphi));
        // }
        if(data.selected_phi[l][i2]){
          htheta->Fill(in_angle(pi-(mtheta)));
          hist_filled=true;
          // }
          // }
        }
      }
      
                
      if(hist_filled==false){
        data.selected_phi[l][i2] = false;
      }
    }
  }
  
}

Double_t DATA::convert_north(Double_t mphi){
  mphi = mphi-phi_rotation_north_rad;
  if(mphi<0){
    mphi = 2*pi+mphi;
  }
  return mphi;
}

void DATA::read_tree(){
  f = new TFile(file_name.c_str(),"READ");
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


 }




