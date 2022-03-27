#define memorydepth 10000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void create_fft(Int_t user_value = -1){

  Int_t chn = 3;

  TFile* f;
  TTree* t1;
  ADC_DATA ch;
  TBranch * bch;
  TCanvas* c;
  TF1* fga;

  Int_t n_cut = 7500;
  Double_t maxFreq = 1.e-6/(4e-9); // 1.e-6 to MHz, maxFrequency = npts / (time range = npts*4ns); 

  Int_t ntotal = 0;

 


  f = new TFile("analyzed.root","READ");
  t1 = (TTree*)f->Get("t1");
  bch = t1->GetBranch(Form("Ch%d",chn));
  Int_t crefXara = 2;
  if(chn == 1) crefXara = 2; // Ch1 (wave0) is C2XArapuca
  if(chn == 2) crefXara = 1; // Ch1 (wave0) is C2XArapuca
  string channel = Form("C%dMiniArapucaA4ch1",crefXara);
  bch->SetAddress(&ch);
  TH1D* hcut = new TH1D("hcut","hcut",n_cut,0,n_cut*4);
  TH1D* h = new TH1D("h","h",memorydepth,0,memorydepth*4);
  TH1D *cFFT_cut = new TH1D(Form("%s 30#mus",channel.c_str()),Form("%s 30#mus",channel.c_str()),n_cut,0,maxFreq) ;
  TH1D *cFFT = new TH1D(Form("%s 40#mus",channel.c_str()),Form("%s 40#mus",channel.c_str()),memorydepth,0,maxFreq) ;
  THStack *hs = new THStack(Form("hstack%d",crefXara),Form("hstack%d",crefXara));
  Int_t start_for;
  Int_t end_for;
  Int_t nevents=0;
  Int_t nx = 0;
  if(user_value == -1){
    start_for = 0;
    end_for = t1->GetEntries();
    // nevents = end_for;
    end_for = 20000; //fast and furious
    nevents = end_for;
  }
  else{
    start_for = user_value;
    end_for = user_value+1;
    nevents=1;
  }
  // for(Int_t j = 0; j<1; j++){
  for(Int_t j = start_for; j<end_for; j++){
    bch->GetEvent(j);
    TH1* hfft=0;
    for(Int_t k = 0; k<n_cut; k++){
      hcut->SetBinContent(k+1,ch.wvf[k]);
    }
    for(Int_t k = 0; k<memorydepth; k++){
      h->SetBinContent(k+1,ch.wvf[k]);
    }

    nx = memorydepth;
    TVirtualFFT *fft = TVirtualFFT::FFT(1,&nx,"R2C P");
    TVirtualFFT::SetTransform(fft);
    hfft = h->FFT(hfft, "MAG");
    for(Int_t k = 0; k<memorydepth; k++){
      cFFT->AddBinContent(k+1,hfft->GetBinContent(k));

    }

    delete hfft;
    hfft=0;

    nx = n_cut;
    TVirtualFFT *fft_cut = TVirtualFFT::FFT(1,&nx,"R2C P");
    TVirtualFFT::SetTransform(fft_cut);
    hfft = hcut->FFT(hfft, "MAG");
    for(Int_t k = 0; k<memorydepth; k++){
      cFFT_cut->AddBinContent(k+1,hfft->GetBinContent(k));

    }

    delete hfft;
    hfft=0;

    
  }

  cFFT->SetEntries(nevents);
  cFFT_cut->SetEntries(nevents);
  // cFFT->SetStats(kFALSE);

  if(crefXara==1){
    cFFT_cut->SetLineColor(kRed);
    cFFT->SetLineColor(kBlue);
  }
  if(crefXara==2){
    cFFT_cut->SetLineColor(kViolet);
    cFFT->SetLineColor(kBlack);
  }
  // c = new TCanvas("c1","c1",1920,0,1920,1080);
  // c->SetLogy(1);
  // cFFT->Draw();
  // cFFT_cut->Draw("SAME");

  // cFFT->GetYaxis()->SetTitle("Magnitude");
  // cFFT->GetXaxis()->SetTitle("Frequency (MHz)");
  
  
  TCanvas *cstack = new TCanvas(Form("cstack%d",crefXara),Form("cstack%d",crefXara),1920,0,1920,1080);
  cstack->SetLogy(1);
  // cFFT->Scale(1./cFFT->Integral("width"));
  // cFFT_cut->Scale(1./cFFT_cut->Integral("width"));
  hs->Add(cFFT);
  hs->Add(cFFT_cut);
  hs->Draw("nostack");
  hs->GetYaxis()->SetTitle("Magnitude");
  hs->GetXaxis()->SetTitle("Frequency (MHz)");
  hs->GetXaxis()->SetRangeUser(0,maxFreq/2.);
  cstack->BuildLegend();
  TFile *fout = new TFile("graphs/ffts_no_filter.root","UPDATE");
  fout->WriteObject(hs,hs->GetName(),"TObject::kOverwrite");
  

  
}
