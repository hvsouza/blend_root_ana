

void plot_together(){
  TFile *fmu = new TFile("averaged_waveforms.root","READ");
  TFile *fspe = new TFile("sphe_waveform_Ch1.root","READ");
  TGraph *gmu = (TGraph*)fmu->Get("average_normalized_ch1");
  TGraph *gspe = (TGraph*)fspe->Get("mean_ch1");
  Int_t n = gmu->GetN();
  Double_t dtime = 4;
  TH1D *hmu = new TH1D("hmu","hmu",n,0,n*dtime);
  TH1D *hspe = new TH1D("hspe","hspe",n,0,n*dtime);

  Double_t normsphe = -1e12
  for(Int_t i = 0; i<n; i++){
    hmu->SetBinContent(i+1,*(gmu->GetY()+i));
    hspe->SetBinContent(i+1,*(gspe->GetY()+i));
    Double_t tmp = *(gspe->GetY()+i);
    if(tmp>=normsphe){
      normsphe = tmp;
    }
  }

  for(Int_t i = 0; i<n; i++){
    hspe->SetBinContent(i+1,*(gspe->GetY()+i));
  }

  hmu->Draw();

  
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
