#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void diff_signal(){

  Double_t dtime = 4;
  Double_t freq = 250;
  TFile *fsphe = new TFile("sphe_waveforms_Ch1.root","READ");
  TGraph *gsphe = (TGraph*)fsphe->Get("mean_ch1");
  
  TH1D *hsphe = new TH1D("hsphe","hsphe",memorydepth*dtime,0,dtime*memorydepth);
  TH1D *hsignal = new TH1D("hsignal","hsignal",memorydepth*dtime,0,dtime*memorydepth);
  TH1D *hsignal_c = new TH1D("hsignal_c","hsignal_c",memorydepth*dtime,0,dtime*memorydepth);
  TH1D *hres = new TH1D("hres","hres",memorydepth*dtime,0,dtime*memorydepth);

   TF1 *func_fit = new TF1("func_fit","([0]*exp(-(x-[2])/[1])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([1]*[1])))*TMath::Erfc((([2]-x)/[3]+[3]/[1])/TMath::Power(2,0.5))+([4]*exp(-(x-[2])/[5])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([5]*[5])))*TMath::Erfc((([2]-x)/[3]+[3]/[5])/TMath::Power(2,0.5))+([6]*exp(-(x-[2])/[7])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([7]*[7])))*TMath::Erfc((([2]-x)/[3]+[3]/[7])/TMath::Power(2,0.5)) + [8]*TMath::Erfc((([2]-x)/[3])/TMath::Power(2,0.5))",350,20000);
  func_fit->SetNpx(20000);
  func_fit->SetParameters(176,909,324,45.7,-156,1027,25.77,88.9,0);
  func_fit->FixParameter(8,0);
  gsphe->Fit("func_fit","R0");
   

  
  for(Int_t i=0; i<memorydepth*4; i++){
    hsphe->SetBinContent(i+1,func_fit->Eval(i));
  }

  Int_t nphotons = 10;
  
  Double_t stddev = 0.25;

  gRandom->SetSeed(1);
  vector<Double_t> ph_time(nphotons);
  vector<TLine*> ph_line(nphotons);
  Double_t photons_std=20; // std in ns of distribuction
  Double_t signal_trigger = 4000;
  Double_t start_sphe = 0;
  for(Int_t i = 0; i<nphotons; i++){
    ph_time[i] = gRandom->Gaus(signal_trigger,photons_std);
    ph_line[i] = new TLine(ph_time[i], 0., ph_time[i], 160);
    ph_line[i]->SetLineColor(kBlack);
    ph_line[i]->SetLineWidth(2);
    
  }

  Double_t val=0;
  Double_t valnoise=0;
  Int_t shift = 200; // in ns
  for(Int_t i = 0; i<dtime*memorydepth; i++){
    for(Int_t k = 0; k<nphotons; k++){
      Double_t time_of_spe = i - ph_time[k]; // If the event happened at 5000 ns, this is the "zero" of the s.p.e
      if(time_of_spe>=start_sphe){
        val = func_fit->Eval(time_of_spe);
        hsignal->AddBinContent(i+1,val);
        valnoise = gRandom->Gaus(0,stddev);
        hsignal->AddBinContent(i+1,valnoise);
        valnoise = gRandom->Gaus(0,stddev);
        hsignal_c->AddBinContent(i+1,valnoise);

        if(i+1+shift<memorydepth*dtime) hsignal_c->AddBinContent(i+1+shift,-val);
      }
      else{
        valnoise = gRandom->Gaus(0,stddev);
        hsignal->AddBinContent(i+1,valnoise);
        valnoise = gRandom->Gaus(0,stddev);
        hsignal_c->AddBinContent(i+1,valnoise);

      }
      
    }
  
  }

  for(Int_t i = 0; i<dtime*memorydepth; i++){
    hres->SetBinContent(i+1,hsignal->GetBinContent(i+1)-hsignal_c->GetBinContent(i+1));
  }
  

  

  TCanvas *c1 = new TCanvas("c1","testing plataform");
  hsignal->Draw();
  hsignal_c->Draw("same");
  // for (const auto &line : ph_line) line->Draw("same");


  TCanvas *csphe = new TCanvas("csphe","sphe fit");
  gsphe->Draw("ALP");
  func_fit->SetRange(0,dtime*memorydepth);
  func_fit->Draw("SAME");

  TCanvas *cres = new TCanvas("cres","cres");
  hres->Draw();
  
  
}
 


