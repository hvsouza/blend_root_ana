#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void integrate_cut(){

  ANALYZER z("z");
  z.setAnalyzer();
  z.kch = 1;
  int kch = z.kch;

  Double_t spe =  1547.06;
  Double_t spe2 = 3173.21;
  Double_t gain = spe2-spe;
  Double_t deltaplus = 1.3;
  Double_t deltaminus = 1.5;

  //use -30 and selection for 5.8
  // z.getSelection("Ch4.selection == 0");
  z.selectByAmplitude(16,10200,10800,-10,"lower");

  TH1D *h = new TH1D("h","h",50000,0,0);
  ADC_DATA chout;
  chout.Set_npts(memorydepth);
  TFile *fout = new TFile("spe_cut.root","RECREATE");
  TFile *fwvf = new TFile("spe_wvf_cut.root","RECREATE");
  TTree *tout = new TTree("tout","tout");
  tout->Branch("Ch4.", &chout);

  vector<Double_t> avwvf(memorydepth,0);
  vector<Double_t> avtemp(memorydepth,0);
  Int_t navg = 0;
  for(Int_t i = 0; i < z.lev->GetN(); i++){
    int ev = z.lev->GetEntry(i);
    z.getWaveform(ev,1);
    z.applyDenoise(16);
    z.integrate(10320,10440);
    h->Fill(z.temp_charge);
    chout.event = ev;
    chout.selection = 0;
    chout.charge = z.temp_charge;
    z.getWaveform(ev,1);
    if(z.temp_charge < gain*deltaplus && z.temp_charge > gain/deltaminus){
      chout.selection = 1;
      for(Int_t j = 0; j < memorydepth; j++){
        avwvf[j] += z.ch[kch]->wvf[j];
        chout.wvf[j] = z.ch[kch]->wvf[j];
      }
      navg++;
    }
    else{
      for(Int_t j = 0; j < memorydepth; j++){
        chout.wvf[j] = z.ch[kch]->wvf[j];
      }
    }
    tout->Fill();
  }
  for(Int_t j = 0; j < memorydepth; j++){
    avwvf[j] = avwvf[j]/navg;
    avtemp[j] = j*4;
  }
  fout->WriteTObject(h,"analyzed_4");

  TGraph *g = new TGraph(memorydepth, &avtemp[0], &avwvf[0]);
  fwvf->WriteTObject(tout,"t1","TObject::kOverwrite");
  fwvf->WriteTObject(g,"mean_ch4");


}
