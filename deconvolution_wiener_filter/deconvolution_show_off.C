#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void deconvolution_show_off(){

  Double_t dtime = 4;
  // Double_t dtime = 16;
  Double_t freq = 250;
  // Double_t freq = 62.5;
  TFile *fsphe = new TFile("sphe_waveforms_Ch1.root","READ");
  // TFile *fsphe = new TFile("spe_template_62_5MHz.root","READ");
  TGraph *gsphe = (TGraph*)fsphe->Get("mean_ch1");
  // TGraph *gsphe = (TGraph*)fsphe->Get("spe_template_62_5MHz");
  
  TH1D *hsphe = new TH1D("hsphe","hsphe",memorydepth,0,dtime*memorydepth);
  TH1D *hsignal = new TH1D("hsignal","hsignal",memorydepth,0,dtime*memorydepth);
  TH1D *havg = new TH1D("havg","havg",memorydepth,0,dtime*memorydepth);
  TH1D *hres = new TH1D("hres","hres",memorydepth,0,dtime*memorydepth);
  TH1D *hexp = new TH1D("hexp","hexp",memorydepth,0,dtime*memorydepth);

  WIENER wpe("wpe",dtime,freq,1e-9,1e6,memorydepth);
  WIENER wsignal("wsignal",dtime,freq,1e-9,1e6,memorydepth);
  WIENER wdec("wdec",dtime,freq,1e-9,1e6,memorydepth);
  WIENER wnoise("wnoise",dtime,freq,1e-9,1e6,memorydepth);
  WIENER wfilter("wfilter",dtime,freq,1e-9,1e6,memorydepth);

  wdec.baseline = 4400;
  Bool_t startfilling = false; // to avoid the first signal oscillation
  Int_t start_sphe = 0;
  // Double_t sphe_pretrigger_ticks = 0;
  Double_t sphe_pretrigger_ticks = 68;
  // Double_t sphe_pretrigger_ticks = 62;
  Double_t offset = 5800-sphe_pretrigger_ticks*dtime;
  Double_t sphe_final_time = 6200;
  // Double_t sphe_final_time = 625*dtime;


  wsignal.maxBin = offset/dtime+3;
  Double_t sphe_integral = 0;
  for(Int_t i=0; i<memorydepth; i++){
    Double_t sphe_val = *(gsphe->GetY()+i);
    if(i>sphe_pretrigger_ticks && sphe_val>=0 && startfilling==false){
      startfilling = true;
      start_sphe = i;
    }
    if(startfilling){
      if(i<sphe_final_time/dtime){        
        hsphe->SetBinContent(i+1,sphe_val);
        // cout << sphe_val << endl;
        if(sphe_val>0)
          sphe_integral+=sphe_val*dtime;
      }
    }
    // hsphe->SetBinContent(i+1,func_fit->Eval(i*dtime));
  }

  Double_t stddev = 0.25;
  wnoise.build_noise(stddev,0);
  wpe.fft(hsphe);
  wfilter.wienerGenFilter(wpe,wnoise,wsignal);

  gRandom->SetSeed(1);
  vector<Double_t> ph_time;
  vector<TLine *> ph_line;
  Double_t fast = 7, slow = 1400;
  Double_t prob_fast = 0.23;
  Double_t prob_slow = 1-prob_fast;
  Int_t nphotons = 1000;
  Double_t signal_integral = 0;

  Double_t time_of_change = 0;

  TF1 *fplar = new TF1("fplar","[0]*exp(-x/[1])/[1]+[2]*exp(-x/[3])/[3]",0,10000);
  fplar->SetParameters(prob_fast,fast,prob_slow,slow);

  vector<Double_t> ratio_det;
  vector<Double_t> ratio_int;
  vector<Double_t> ratio_count;
  Double_t photons_count = 0;
  Double_t func_int = 0;
  for(Int_t event = 0; event<1; event++){
    hsignal->Reset();
    // hexp->Reset();
    // havg->Reset();
    // hres->Reset();
    // nphotons = gRandom->Poisson(200);
    // nphotons = gRandom->Uniform(200,1000);
    nphotons = 500;
    ph_time.resize(nphotons,0);
    ph_line.resize(nphotons,0);
      
    for(Int_t i = 0; i<nphotons; i++){
      Double_t dice1 = gRandom->Uniform(0,1);
      if(dice1<prob_fast){
        ph_time[i] = gRandom->Exp(fast)+offset;
      }
      else{
        ph_time[i] = gRandom->Exp(slow)+offset;
      }

      ph_line[i] = new TLine(ph_time[i], 0., ph_time[i], 160);
      ph_line[i]->SetLineColor(kBlack);
      ph_line[i]->SetLineWidth(2);
    
    }


    Double_t val=0;
    Double_t valnoise=0;
    for(Int_t i = 0; i<memorydepth; i++){
      for(Int_t k = 0; k<nphotons; k++){
        Double_t time_of_spe = i*dtime + sphe_pretrigger_ticks*dtime- ph_time[k]; // If the event happened at 5000 ns, this is the "zero" of the s.p.e
        hexp->AddBinContent(i+1,419*TMath::Gaus(i*dtime,ph_time[k],20.5,false));
        // hexp->AddBinContent(i+1,138*TMath::Gaus(i*dtime,ph_time[k],45.7,false));
        if(time_of_spe<sphe_final_time && time_of_spe>=start_sphe*dtime){
          val = gsphe->Eval(time_of_spe);
          // val = gsphe->Eval(time_of_spe/1.e3);
          hsignal->AddBinContent(i+1,val);
          havg->AddBinContent(i+1,val);
        }
        else continue;
      
      }
      valnoise = gRandom->Gaus(0,stddev);
      hsignal->AddBinContent(i+1,valnoise);
      havg->AddBinContent(i+1,valnoise);
      hexp->AddBinContent(i+1,valnoise);

    }
    wsignal.fft(hsignal);
    wdec.frequency_deconv(wsignal,wfilter);
    
    Bool_t check_change = true;
    signal_integral = 0;
    for(Int_t i=0; i<memorydepth; i++){
      Double_t signal_val = hsignal->GetBinContent(i+1);
      hres->AddBinContent(i+1,wdec.hwvf->GetBinContent(i+1));
      if(signal_val>0)signal_integral+=signal_val*dtime;
      if(i>6000/4 && check_change){
        if(signal_val<=0){
          time_of_change = i*dtime;
          check_change = false;
        }
      }
    }
    cout << event << endl;
    photons_count = 0;
    for(Int_t i=0; i<nphotons; i++){
      if(ph_time[i]<time_of_change){
        photons_count+=1;
      }
    }
    func_int = fplar->Integral(0,time_of_change-offset);
    ratio_count.push_back(photons_count/nphotons);
    ratio_int.push_back(func_int);
    ratio_det.push_back(signal_integral/sphe_integral/nphotons);
 
  }

  

  cout << "sphe: " << sphe_integral << endl;
  cout << "signal: " << signal_integral << endl;
  cout << "n: " << signal_integral/sphe_integral << endl;
  cout << "start: " << offset << endl;
  cout << "finish: " << time_of_change << endl;
  cout << "integrated: " << func_int << endl;
  cout << "count photons: " << photons_count << endl;
  
  Double_t maxVal = hsignal->GetMaximum();
  // hsignal->Scale(1/maxVal);


  // wdec.shift_waveform(hres,offset/4);
  maxVal = hsphe->GetMaximum();
  // hsphe->Scale(1/maxVal);

  TCanvas *c1 = new TCanvas("c1","testing plataform");
  hsignal->Draw();
  // fplar->Draw();
  TCanvas *cexp = new TCanvas("cexp","expected waveform");
  hexp->Draw();
  for (const auto &line : ph_line) line->Draw("same");


  //  Double_t baseline = 0;
  // Double_t auxbaseline = 0;
  // for(Int_t i=0; i<5000/dtime; i++){
  //   baseline+=hres->GetBinContent(i+1);
  //   auxbaseline+=1;
  // }
  // baseline=baseline/auxbaseline;
  // for(Int_t i=0; i<memorydepth; i++){
  //   hres->SetBinContent(i+1,hres->GetBinContent(i+1)-baseline);
  // }
 
  
  TCanvas *c2 = new TCanvas("c2","processed");
  
  hres->Draw();
  hexp->SetLineColor(kRed);
  hexp->Draw("SAME");


  TCanvas *c3 = new TCanvas("c3","ratio by integral");
  TGraph *gr1 = new TGraph(ratio_int.size(),&ratio_count[0], &ratio_int[0]);
  gr1->Draw("AP");
  TCanvas *c4 = new TCanvas("c4","ratio by detection");
  TGraph *gr2 = new TGraph(ratio_int.size(),&ratio_count[0], &ratio_det[0]);
  gr2->Draw("AP");

}
 


