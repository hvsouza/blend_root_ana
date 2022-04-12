#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void create_histograms(){
  vector<Double_t> volts = {5,7.5,10,12.5,15,17.5,20,22,25,27.5,30};
  Int_t n = volts.size();
//   vector<string> filesC2 = {"C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl7p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl10V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl12p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl15V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl17p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl20V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl22p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl25V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl27p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl30V00000"};
  vector<string> filesC2 = {"C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl7p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl10V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl12p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl15V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl17p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl20V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl22p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl25V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl27p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl30V00000"};

  string conc = "/analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);
  vector<TH1D*> h(n);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  TCanvas *c1 = new TCanvas();
  vector<TF1*> fga(n);
  Double_t charge = 0;
  c1->Divide(4,3);
  c1->Update();

  Double_t sphe = 0.3;  // V.ns
  
  for(Int_t i = 0; i<n; i++){
    filesC2[i] = filesC2[i]+conc;
    f[i] = new TFile(filesC2[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch("Ch1");
    bch[i]->SetAddress(&ch[i]);
    h[i] = new TH1D(Form("h%d",i),Form("h%d",i),200,-2/sphe,3/sphe);
    // h[i] = new TH1D(Form("h%d",i),Form("h%d",i),200,-0.01,0.01);
    
    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
      bch[i]->GetEvent(j);
      charge = 0;
      // for(Int_t k = 0/4; k<1000/4; k++){
      //  charge = ch[i].wvf[k];
      //  h[i]->Fill(charge);
      // }
      for(Int_t k = 1730/4; k<2000/4.; k++){
       charge += ch[i].wvf[k];
      }
      h[i]->Fill(charge/sphe);
    }
    c1->cd(i+1);
    h[i]->Draw();
//     h[i]->GetXaxis()->SetRangeUser(-0.04,0.04);
    string hist_name = Form("C2XArapuca A1ch1 - LED %.1f V",volts[i]);
    h[i]->SetNameTitle(hist_name.c_str(),hist_name.c_str());
    h[i]->GetYaxis()->SetTitle("# of events");
    h[i]->GetXaxis()->SetTitle("Photo-electrons");
    fga[i] = new TF1(Form("f%d",i),"gaus(0)",-1,1);
    fga[i]->SetParameters(0.1,0.1,0.1);
    h[i]->Fit(fga[i]);
    fga[i]->Draw("SAME");
  }
  
  vector<Double_t> avgs(n);
  vector<Double_t> sigma_mu(n);
  vector<Double_t> Ersigma_mu(n);
  vector<Double_t> Eravgs(n);
  vector<Double_t> ErSTD(n);
  for(Int_t i = 0; i<n; i++){
    avgs[i] = fga[i]->GetParameter(1);
    ErSTD[i] = fga[i]->GetParameter(2);
    Eravgs[i] = sqrt(pow(fga[i]->GetParError(1),2)+pow(ErSTD[i],2)/h[i]->GetEntries());

    sigma_mu[i] = 100*ErSTD[i]/avgs[i];
    Ersigma_mu[i] = sigma_mu[i]*sqrt(pow(fga[i]->GetParError(2)/ErSTD[i],2)+pow(fga[i]->GetParError(1)/avgs[i],2));
  }
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas *c2 = new TCanvas();

  TGraphErrors *gmeans = new TGraphErrors(n,&volts[0],&avgs[0],0,&Eravgs[0]);
//   TGraphErrors *gSTDs = new TGraphErrors(n,&volts[0],&avgs[0],0,&ErSTD[0]);
  gmeans->GetXaxis()->SetTitle("LED Voltage (V)");
  gmeans->GetYaxis()->SetTitle("Average photo-electrons");
  gmeans->SetMarkerStyle(21);
  gmeans->SetMarkerSize(0.5);
//   gSTDs->SetMarkerSize(0);
//   gSTDs->SetLineWidth(2);
//   gSTDs->Draw("AP []");
  gmeans->Draw("AP");
  gmeans->Fit("pol1");
  
  TCanvas *c3 = new TCanvas();

  TGraphErrors *gstdmu = new TGraphErrors(n,&avgs[0],&sigma_mu[0],&Eravgs[0],&Ersigma_mu[0]);
  TF1 *myfunc = new TF1("myfunc","sqrt([0]*[0]+[1]*[1]/x+pow([2]/x,2))",0,35);
  myfunc->SetParameters(1,2,3);
  gstdmu->GetXaxis()->SetTitle("#mu (photo-electrons)");
  gstdmu->GetYaxis()->SetTitle("#sigma/#mu (%)");
  gstdmu->SetMarkerStyle(21);
  gstdmu->SetMarkerSize(0.5);
  gstdmu->Draw("AP");
  gstdmu->Fit("myfunc");

  


  
  
  
  
}
