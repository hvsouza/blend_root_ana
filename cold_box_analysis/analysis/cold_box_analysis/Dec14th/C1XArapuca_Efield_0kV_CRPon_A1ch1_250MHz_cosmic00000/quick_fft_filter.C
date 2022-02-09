#define memorydepth 2502
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void quick_fft_filter(int user_value = 0){
  Int_t n = memorydepth;
  
  TFile *f1 = new TFile("analyzed.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  ADC_DATA ch;
  TBranch *bch = t1->GetBranch("Ch1");
  bch->SetAddress(&ch);
  Int_t nentries = t1->GetEntries();
  
  
  TH1D *hsignal = new TH1D("hsignal", "hsignal", memorydepth, 0, memorydepth);   
  TH1D *hfiltered = new TH1D("hfiltered", "hfiltered", memorydepth, 0, memorydepth);   
  
  t1->GetEntry(user_value);
  
  //Fill the histogram with function values
  for (Int_t i=0; i<memorydepth; i++){
    hsignal->SetBinContent(i+1, ch.wvf[i]);
  }
  
  TCanvas *c1 = new TCanvas();
  
  TH1 *hm =0;
  TH1 *hcut =0;
  TH1 *hcutre =0;
  TH1 *hcutim =0;
  
  TVirtualFFT::SetTransform(0);
  hm = hsignal->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  
  
  TH1 *holdmag = (TH1*)hm->Clone("holdmag");
  holdmag->Draw();
  TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();
  
  //Use the following method to get the full output:
  Double_t *re_full2 = new Double_t[n];
  Double_t *im_full2 = new Double_t[n];
  fft2->GetPointsComplex(re_full2,im_full2);
  
  
  
  
  Double_t *re_final = new Double_t[n];
  Double_t *im_final = new Double_t[n];
  TF1 *filter = new TF1("filter","TMath::Gaus(x,[0],[1])",0,n);	// A gaussian filter
  filter->SetParameters(0,55); // USE THIS!!! -> center = 0, cut = 180
  
  for(Int_t i = 0; i<n; i++){
    
    re_final[i] = filter->Eval(i)*re_full2[i];
    im_final[i] = filter->Eval(i)*im_full2[i];
  }
  
  
  
  //Now let's make a backward transform:
  TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_final->SetPointsComplex(re_final,im_final);
  fft_final->Transform();
  TH1 *hfinal = 0;
  
  
  
  //Let's look at the output
  hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
  hfinal->Scale(1./memorydepth);
  
  TVirtualFFT::SetTransform(0);
  hcut = hfinal->FFT(hcut, "MAG");
  hcut->SetTitle("Magnitude of the 1st transform");
  hcut->Draw("SAME");
  hcut->SetLineColor(kRed);
  
  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  hfinal->Draw("hist");
  hsignal->Draw("SAME hist");
  hfinal->SetLineWidth(2);
  
  hfinal->SetLineColor(kRed);
  hfinal->GetXaxis()->SetLabelSize(0.05);
  hfinal->GetYaxis()->SetLabelSize(0.05);
  
  
  
}
