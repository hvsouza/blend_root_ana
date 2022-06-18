#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void normal_deconvolution(){

  TFile *fsphe = new TFile("sphe_waveforms_Ch1_no_filter.root","READ");
  TGraph *gsphe = (TGraph*)fsphe->Get("mean_ch1");

  TFile *favg = new TFile("averaged_waveforms_no_filter.root","READ");
  TGraph *gavg = (TGraph*)favg->Get("average_normalized_ch1");

  TFile *f1 = new TFile("analyzed_no_filter.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");

  ADC_DATA ch1;
  TBranch *bch = t1->GetBranch("Ch1");
  bch->SetAddress(&ch1);

  
  TH1D *hsphe = new TH1D("hsphe","hsphe",memorydepth,0,4*memorydepth);
  TH1D *havg = new TH1D("havg","havg",memorydepth,0,4*memorydepth);
  TH1D *hsignal = new TH1D("hsignal","hsignal",memorydepth,0,4*memorydepth);
  TH1D *hres = new TH1D("hres","hres",memorydepth,0,4*memorydepth);

  WIENER wpe("wpe");
  WIENER wavg("wavg");
  WIENER wsignal("wsignal");
  WIENER wback("wback");
  WIENER wdec("wdec");

  TF1 *func_fit = new TF1("func_fit","([0]*exp(-(x-[2])/[1])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([1]*[1])))*TMath::Erfc((([2]-x)/[3]+[3]/[1])/TMath::Power(2,0.5))+([4]*exp(-(x-[2])/[5])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([5]*[5])))*TMath::Erfc((([2]-x)/[3]+[3]/[5])/TMath::Power(2,0.5))+([6]*exp(-(x-[2])/[7])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([7]*[7])))*TMath::Erfc((([2]-x)/[3]+[3]/[7])/TMath::Power(2,0.5)) + [8]*TMath::Erfc((([2]-x)/[3])/TMath::Power(2,0.5))",350,20000);
  func_fit->SetNpx(20000);
  func_fit->SetParameters(176,909,324,45.7,-156,1027,25.77,88.9,0);
  func_fit->FixParameter(8,0);
  gsphe->Fit("func_fit","R0");
   

  Bool_t startfilling = false; // to avoid the first signal oscillation
  for(Int_t i=0; i<memorydepth; i++){
    if(i*4>200 && *(gsphe->GetY()+i)>=0){
      startfilling = true;
    }
    if(startfilling) hsphe->SetBinContent(i+1,*(gsphe->GetY()+i)); 
    // hsphe->SetBinContent(i+1,func_fit->Eval(i*4));
    havg->SetBinContent(i+1,*(gavg->GetY()+i));
    // hsignal->SetBinContent(i+1,ch1.wvf[i]);
  }

  // wpe.shift_waveform(hsphe,3500);
  // wpe.shift_waveform(hsignal,2500);

  wpe.fft(hsphe);
  wback.backfft(wpe);

  for(Int_t event = 0; event<t1->GetEntries(); event++){
  // for(Int_t event = 0; event<100; event++){
     bch->GetEvent(event);
     // if(ch1.selection==0){
       if(ch1.fprompt<0.5){
         if(ch1.charge>2.2e6 && ch1.charge<12e6 && ch1.peak<5000){
           if(ch1.wvf[12000/4]<-1500) continue; // there is one event with saturated undershoot, I want to avoid it..
           for(Int_t i=0; i<memorydepth; i++){
             hsignal->SetBinContent(i+1,ch1.wvf[i]);
           }
           wsignal.fft(hsignal);
           wdec.deconvolve(wsignal,wpe,2);
           for(Int_t i=0; i<memorydepth; i++){
             hres->AddBinContent(i+1,wdec.hwvf->GetBinContent(i+1));
           }
           
           
         }
       }
     // }
  }

  
  Double_t baseline = 0;
  Double_t auxbaseline = 0;
  for(Int_t i=0; i<5000/4.; i++){
    baseline+=hres->GetBinContent(i+1);
    auxbaseline+=1;
  }
  baseline=baseline/auxbaseline;
  for(Int_t i=0; i<memorydepth; i++){
    hres->SetBinContent(i+1,hres->GetBinContent(i+1)-baseline);
  }
           
  


  TCanvas *cpe = new TCanvas("cpe","sphe");
  // gsphe->Draw("");
  // func_fit->SetRange(0,memorydepth*4);
  // func_fit->Draw("SAME");
  // hsignal->Draw("");
  hsphe->Draw("");

  

  TCanvas *c = new TCanvas("c","ffts");
  wpe.hfft->Draw();
  wavg.hfft->Draw("SAME");

  // wdec.deconvolve(wavg,wpe,2);

  TCanvas *c2 = new TCanvas("c2","deconvolved");
  hres->Draw();


 
  
}
