#define memorydepth 1000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void deconvolution_show_off(){

  Double_t dtime = 16;
  // TFile *fsphe = new TFile("sphe_waveforms_Ch1.root","READ");
  TFile *fsphe = new TFile("spe_template_62_5MHz.root","READ");
  // TGraph *gsphe = (TGraph*)fsphe->Get("mean_ch1");
  TGraph *gsphe = (TGraph*)fsphe->Get("spe_template_62_5MHz");
  
  TH1D *hsphe = new TH1D("hsphe","hsphe",memorydepth,0,dtime*memorydepth);
  TH1D *hsignal = new TH1D("hsignal","hsignal",memorydepth,0,dtime*memorydepth);
  TH1D *havg = new TH1D("havg","havg",memorydepth,0,dtime*memorydepth);
  TH1D *hres = new TH1D("hres","hres",memorydepth,0,dtime*memorydepth);

  WIENER wpe("wpe");
  WIENER wsignal("wsignal");
  WIENER wdec("wdec");
  WIENER wnoise("wnoise");
  WIENER wfilter("wfilter");

  Bool_t startfilling = false; // to avoid the first signal oscillation
  Int_t start_sphe = 0;
  for(Int_t i=0; i<memorydepth; i++){
    // if(i*dtime>200 && *(gsphe->GetY()+i)>=0 && startfilling==false){
    if(i>68 && *(gsphe->GetY()+i)>=0 && startfilling==false){
      startfilling = true;
      start_sphe = i;
    }
    if(startfilling){
      if(i<625)
        hsphe->SetBinContent(i+1,*(gsphe->GetY()+i));
    }
    // hsphe->SetBinContent(i+1,func_fit->Eval(i*dtime));
  }

  Double_t stddev = 0.25;
  wnoise.build_noise(stddev,0);
  wpe.fft(hsphe);
  wfilter.wienerGenFilter(wpe,wnoise,wsignal);

  gRandom->SetSeed(0);
  vector<Double_t> ph_time;
  vector<TLine *> ph_line;
  Double_t fast = 7, slow = 1400;
  Double_t prob_fast = 0.23;
  Double_t prob_slow = 1-prob_fast;
  Double_t offset = 5800-62*dtime
    ;
  // Double_t sphe_final_time = 6300;
  Double_t sphe_final_time = 625*dtime;

  
  for(Int_t event = 0; event<20; event++){
    hsignal->Reset();
    // Int_t nphotons = gRandom->Poisson(200);
    Int_t nphotons = 200;
    ph_time.resize(nphotons);
    ph_line.resize(nphotons);
      
    for(Int_t i = 0; i<nphotons; i++){
      Double_t dice1 = gRandom->Uniform(0,1);
      if(dice1<prob_fast){
        ph_time[i] = gRandom->Exp(fast)+offset;
      }
      else{
        ph_time[i] = gRandom->Exp(slow)+offset;
      }

      ph_line[i] = new TLine(ph_time[i]+start_sphe*dtime, 0., ph_time[i]+start_sphe*dtime, 0.1);
      ph_line[i]->SetLineColor(kBlack);
      ph_line[i]->SetLineWidth(2);
    
    }


    Double_t val=0;
    Double_t valnoise=0;
    for(Int_t i = 0; i<memorydepth; i++){
      for(Int_t k = 0; k<nphotons; k++){
        Double_t time_of_spe = i*dtime - ph_time[k]; // If the event happened at 5000 ns, this is the "zero" of the s.p.e
        if(time_of_spe<sphe_final_time && time_of_spe>=start_sphe*dtime){
          val = gsphe->Eval(time_of_spe/1.e3,0,"S");
          hsignal->AddBinContent(i+1,val);
          havg->AddBinContent(i+1,val);
        }
        else continue;
      
      }
      valnoise = gRandom->Gaus(0,stddev);
      hsignal->AddBinContent(i+1,valnoise);
      havg->AddBinContent(i+1,valnoise);

    }
    wsignal.fft(hsignal);
    wdec.frequency_deconv(wsignal,wfilter);

    for(Int_t i=0; i<memorydepth; i++){
      hres->AddBinContent(i+1,wdec.hwvf->GetBinContent(i+1));
    }
    cout << event << endl;
  }

  Double_t maxVal = hsignal->GetMaximum();
  // hsignal->Scale(1/maxVal);


  // wpe.shift_waveform(hsphe,6000/4.);
  maxVal = hsphe->GetMaximum();
  // hsphe->Scale(1/maxVal);

  TCanvas *c1 = new TCanvas("c1","testing plataform");
  hsignal->Draw();
  for (const auto &line : ph_line) line->Draw("same");


   Double_t baseline = 0;
  Double_t auxbaseline = 0;
  for(Int_t i=0; i<5000/dtime; i++){
    baseline+=hres->GetBinContent(i+1);
    auxbaseline+=1;
  }
  baseline=baseline/auxbaseline;
  for(Int_t i=0; i<memorydepth; i++){
    hres->SetBinContent(i+1,hres->GetBinContent(i+1)-baseline);
  }
 
  
  TCanvas *c2 = new TCanvas("c2","processed");
  
  hres->Draw();

  
  

}
 


