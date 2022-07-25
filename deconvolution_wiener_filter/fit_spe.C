#define memorydepth 5000


void fit_spe(){

  TFile *fsphe = new TFile("sphe_waveforms_Ch1_no_filter.root","READ");
  TGraph *gsphe = (TGraph*)fsphe->Get("mean_ch1");
  TH1D *hsphe = new TH1D("hsphe","hsphe",memorydepth,0,4*memorydepth);

  TF1 *func_fit = new TF1("func_fit","([0]*exp(-(x-[2])/[1])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([1]*[1])))*TMath::Erfc((([2]-x)/[3]+[3]/[1])/TMath::Power(2,0.5))+([4]*exp(-(x-[2])/[5])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([5]*[5])))*TMath::Erfc((([2]-x)/[3]+[3]/[5])/TMath::Power(2,0.5))+([6]*exp(-(x-[2])/[7])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([7]*[7])))*TMath::Erfc((([2]-x)/[3]+[3]/[7])/TMath::Power(2,0.5)) + [8]*TMath::Erfc((([2]-x)/[3])/TMath::Power(2,0.5))",350,20000);
  func_fit->SetNpx(20000);
  func_fit->SetParameters(176,909,324,45.7,-156,1027,25.77,88.9,0);
  func_fit->FixParameter(8,0);
  gsphe->Fit("func_fit","R0");
   
  for(Int_t i=0; i<memorydepth; i++){
    hsphe->SetBinContent(i+1,*(gsphe->GetY()+i));
  }


  TCanvas *cpe = new TCanvas("cpe","sphe");
  gsphe->Draw("");
  func_fit->SetRange(0,memorydepth*4);
  func_fit->Draw("SAME");
  
 
  
}
