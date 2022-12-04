#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void showPoisson(){
  TFile *f = new TFile("sphe_histograms_italy.root","READ");
  Calibration Cal;

  Cal.dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)

  Cal.rebin = 1;
  Cal.make_free_stddevs = true;
    Cal.rootFile = "SPE.root";
    Cal.searchParameters("analyzed_1",2,true);

  TH1D *h = (TH1D*)gDirectory->Get(Cal.histogram_name.c_str());
  h->Rebin(Cal.rebin);
  h->Scale(1./h->Integral());

  Double_t total_events = h->Integral("width");
  Double_t zero_events = Cal.faux->Integral(Cal.xmin, Cal.xmax);
  Double_t prob_zeros = zero_events/total_events;
  // cout << total_events << " " << zero_events << " " << prob_zeros << endl;
  //

  Double_t width = h->GetBinWidth(1);
  Double_t startbin = h->GetBinCenter(1);

  // bin_center = (nbins-1)*width + startbin;
  Int_t npois = 3;
  Int_t bin_n_peak = round(npois*Cal.mean1 - startbin)/width + 1;

  Double_t lambda1 = -TMath::Log(prob_zeros);
  Double_t lambda2 = h->GetMean()/Cal.mean1;
  Double_t lambda = lambda1;


  Double_t peak_by_formula1 = Cal.startpoint/(pow(Cal.poisson(lambda1,npois-1)/Cal.poisson(lambda1,npois),npois-2));
  Double_t peak_by_formula2 = Cal.startpoint/(pow(Cal.poisson(lambda2,npois-1)/Cal.poisson(lambda2,npois),npois-2));
  Double_t diff1 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula1);
  Double_t diff2 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula2);

  // cout << h->GetBinContent(bin_n_peak) << " " << peak_by_formula1 << " " << peak_by_formula2 <<  " " << diff1 << " " << diff2 << endl;
  if(diff1 <= diff2) lambda = lambda1;
  else lambda = lambda2;
  // lambda = lambda1;

  Double_t poisson2 = Cal.poisson(lambda,2);
  Double_t poisson3 = Cal.poisson(lambda,3);
  Double_t poisson_ratio = poisson2/poisson3;

  Double_t totalEntries = h->GetEntries();
  Double_t events0 = prob_zeros*totalEntries;
  Int_t np = 8;
  vector<Double_t> eventsn(np);
  vector<Double_t> poissonn(np);
  for (Int_t i = 0; i < np; i++) {
    poissonn[i] = Cal.poisson(lambda,i);
    eventsn[i] = poissonn[i]*totalEntries;
  }

  vector<TF1*> fn(np);

  Cal.faux->SetRange(Cal.xmin,Cal.xmax);
  cout << lambda << " " << lambda1 << " " << lambda2 << poisson_ratio << endl;

  for (Int_t i = 0; i < np; i++) {
    fn[i] = (TF1*)Cal.faux->Clone(Form("f%d",i+1));
    Double_t prev = prob_zeros;
    if(i == 0){
      fn[i]->SetParameter(0,Cal.peak0);
    }
    else if(i == 1){
      fn[i]->SetParameter(0,Cal.peak1);
      fn[i]->SetParameter(1,Cal.mean1);
    }
    else{
      Double_t ratio = i > 2 ? poisson_ratio : 1;
      Cal.startpoint = Cal.startpoint/ratio;
      fn[i]->SetParameter(0,Cal.startpoint);
      fn[i]->SetParameter(1,Cal.mean1*i);
    }
    fn[i]->SetLineColor(kBlack);
    fn[i]->SetLineWidth(3);
    fn[i]->SetNpx(10000);
    fn[i]->Draw("SAME");

  }





}
