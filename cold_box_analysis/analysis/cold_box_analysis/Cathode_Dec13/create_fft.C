#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void create_fft(Int_t user_value = 0){
  Int_t n = 1;
  vector<Double_t> volts = {0,5,10};
  vector<string> filesC2 = {"C1XArapuca_Efield_0kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_5kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_10kV_b_A4ch2_250MHz_cosmic00000"};
//   vector<string> filesC2 = {"C2XArapuca_Efield_0kV_A1ch1_250MHz_cosmic00000","C2XArapuca_Efield_5kV_A1ch1_250MHz_cosmic00000","C2XArapuca_Efield_10kV_b_A1ch1_250MHz_cosmic00000"};
  string conc = "/analyzed.root";
  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);
  vector<TH1D*> h(n);
  vector<TH1*> hfft(n);
  vector<TCanvas*> c(n);
  vector<TF1*> fga(n);
  
  Double_t maxFreq = 1252.e-6/(5008e-9);
  TH1D *finalFFT = new TH1D("finalFFT","finalFFT",1252,0,maxFreq) ; // times 10^-3 to give in MHz
  vector<TH1D *> cFFT(n);
  Int_t ntotal = 0;
  
  for(Int_t i = 0; i<n; i++){
    filesC2[i] = filesC2[i]+conc;
    f[i] = new TFile(filesC2[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch("Ch1");
    bch[i]->SetAddress(&ch[i]);
    h[i] = new TH1D(Form("h%d",i),Form("h%d",i),1252,0,5008);
    cFFT[i] = new TH1D(Form("hcfft%d",i),Form("hcfft%d",i),1252,0,maxFreq) ;
//     for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
    for(Int_t j = user_value; j<user_value+1; j++){
      bch[i]->GetEvent(j);
    
      for(Int_t k = 0; k<memorydepth; k++){
        h[i]->SetBinContent(k+1,ch[i].wvf[k]);
      }
      
      TVirtualFFT::SetTransform(0);
      hfft[i] = h[i]->FFT(hfft[i], "MAG");
      for(Int_t k = 0; k<memorydepth; k++){
        cFFT[i]->AddBinContent(k+1,hfft[i]->GetBinContent(k));
        finalFFT->AddBinContent(k+1,hfft[i]->GetBinContent(k));
      }
//       cout << hfft[i]->GetEntries() << endl;
    hfft[i]->Reset();
    }
    
    cFFT[i]->SetEntries(t1[i]->GetEntries());
    ntotal+=t1[i]->GetEntries();
    
    c[i] = new TCanvas(Form("c%d",i));
    c[i]->SetLogy(1);
    cFFT[i]->Draw();

    cFFT[i]->SetNameTitle(Form("X-ARAPUCA C2 - %0.f kV",volts[i]),Form("X-ARAPUCA C2 - %0.f kV",volts[i]));
    cFFT[i]->GetYaxis()->SetTitle("Magnitude");
    cFFT[i]->GetXaxis()->SetTitle("Frequency (MHz)");
    
  }
  

  
  TCanvas *cf = new TCanvas("cf","cf");
  cf->SetLogy(1);
  finalFFT->SetEntries(ntotal);
  finalFFT->Draw();
  
  finalFFT->SetNameTitle("Signals FFT C2","Signals FFT C2");
  finalFFT->GetYaxis()->SetTitle("Magnitude");
  finalFFT->GetXaxis()->SetTitle("Frequency (MHz)");
  finalFFT->GetXaxis()->SetRangeUser(0,maxFreq/2.);
  


  
  
  
  
}
