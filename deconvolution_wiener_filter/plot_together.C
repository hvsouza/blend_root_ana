template <class T>
void shift_waveform(T *h, Int_t new_max){
  Int_t npts = h->GetNbinsX();
  Int_t old_max = h->GetMaximumBin();
  Int_t old_ref = old_max - new_max;
  TH1D *htemp = (TH1D*)h->Clone("htemp");
  Double_t temp;
  if(old_ref<0){
    // cout << " case lower" << endl;
    old_ref = npts-(new_max-old_max);
  }
  for(Int_t i = 1; i<npts-(old_ref); i++){
    temp = htemp->GetBinContent(old_ref+i);
    h->SetBinContent(i,temp);
  }
  Int_t aux = 1;
  for(Int_t i = npts-(old_ref); i<=npts; i++){
    temp = htemp->GetBinContent(aux);
    h->SetBinContent(i,temp);
    aux++;
  }
  delete htemp;
}


void plot_together(){
  TFile *fmu = new TFile("averaged_waveforms.root","READ");
  TFile *fspe = new TFile("sphe_waveforms_Ch1.root","READ");
  TGraph *gmu = (TGraph*)fmu->Get("average_normalized_ch1");
  TGraph *gspe = (TGraph*)fspe->Get("mean_ch1");
  Int_t n = gmu->GetN();
  Double_t dtime = 4;
  TH1D *hmu = new TH1D("hmu","hmu",n,0,n*dtime);
  TH1D *hspe = new TH1D("hspe","hspe",n,0,n*dtime);

  Double_t normsphe = -1e12;
  for(Int_t i = 0; i<n; i++){
    hmu->SetBinContent(i+1,*(gmu->GetY()+i));
    Double_t tmp = *(gspe->GetY()+i);
    if(tmp>=normsphe){
      normsphe = tmp;
    }
  }

  for(Int_t i = 0; i<n; i++){
    hspe->SetBinContent(i+1,*(gspe->GetY()+i)/normsphe);
  }

  shift_waveform(hspe,5986/4);

  hspe->SetLineWidth(2);
  hmu->SetLineWidth(2);
  hspe->SetLineColor(kBlue);
  hspe->SetLineColor(kBlack);
  hspe->SetTitle("SPE");
  hmu->SetTitle("Muon");
  hmu->GetXaxis()->SetTitle("Time (ns)");
  hmu->GetYaxis()->SetTitle("Amplitude (A.U.)");
  hmu->Draw();
  hspe->Draw("same");
  
}


/*
  TFile *fmu = new TFile("averaged_waveforms.root","READ");
  TGraph *gmu = (TGraph*)fmu->Get("average_normalized_ch1");
  Int_t n = gmu->GetN();
  Double_t dtime = 4;
  TH1D *hmu = new TH1D("hmu","hmu",n,0,n*dtime);

  for(Int_t i = 0; i<n; i++){
    hmu->SetBinContent(i+1,*(gmu->GetY()+i));
  }

  hmu->Draw();
*/
