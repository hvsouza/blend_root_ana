#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void create_histograms_ch1(){
  string mychannel = "Ch1";
  vector<Double_t> volts = {3,4,4.2,4.4,4.6,4.8,5,6,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30};
  // vector<Double_t> volts = {3,4,4.2,4.4,4.6,4.8,5,6,7.5};
  // vector<Double_t> volts = {10,12.5,15,17.5,20,22.5,25,27.5,30};
  // vector<Double_t> volts = {7.5};
  Int_t n = volts.size();

  vector<string> files(n);
  for(Int_t i=0; i<n; i++){
    if(static_cast<Int_t>(volts[i])-volts[i]!=0) files[i] = Form("%.1fV",volts[i]);
    else files[i] = Form("%.0fV",volts[i]);
  }
  vector<Double_t> ranges = {-4000,4000};

  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);
  vector<TH1D*> h(n);
  vector<TH1D*> hpeak(n);
  TH2D *h2p = new TH2D("h2p","h2p",400,-50,350,1000,0,1.);
  TH2D *h2 = new TH2D("h2","h2",300,0,30,4000,-50,350);
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  TCanvas *cpeak = new TCanvas("cpeak","cpeak",1920,0,1920,1080);
  vector<TF1*> fga(n);
  vector<TF1*> fgap(n);
  Double_t charge = 0;
  Int_t m1=0,m2=0;
  m1=m2=static_cast<Int_t>(sqrt(n));
  // cout << m2 << endl;
  if(m1*m2==n){}
  else if((m2+1)*m1>=n) m2+=1;
  else if((m2+1)*(m1+1)>=n){
    m2+=1;
    m1+=1;
  }
  
  // cout << m1 << " " << m2 << endl;
  c1->Divide(m1,m2);
  c1->Update();

  cpeak->Divide(m1,m2);
  cpeak->Update();


  Double_t sphe = 0.4;  // V.ns
  Double_t one_adc_channel = (2/pow(2.,14));
  sphe = sphe/one_adc_channel; // ADC.ns 
  Double_t increase = 0;
  Int_t nbins_d = 100;
  for(Int_t i = 0; i<n; i++){
    files[i] = files[i]+conc;
    f[i] = new TFile(files[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch(mychannel.c_str());
    bch[i]->SetAddress(&ch[i]);
    hpeak[i] = new TH1D(Form("hpeak%d",i),Form("hpeak%d",i),100,0,0);
    if(i<4){
      ranges = {-5,10};
      nbins = static_cast<Int_t>((ranges[1]-ranges[0])*12);
      h[i] = new TH1D(Form("h%d",i),Form("h%d",i),nbins,0,0);
    }
    else if(i>=4 && i<7){
      ranges = {-5,15};
      nbins = static_cast<Int_t>((ranges[1]-ranges[0])*12);
      h[i] = new TH1D(Form("h%d",i),Form("h%d",i),nbins,0,0);
    }
    else if(i>=7){
      if(i==9) nbins_d+=50;
      if(i==14) nbins_d+=50;
      if(i==16) nbins_d+20;
      ranges = {-10,12+increase};
      nbins = static_cast<Int_t>((ranges[1]-ranges[0])/nbins_d);
      h[i] = new TH1D(Form("h%d",i),Form("h%d",i),nbins_d,0,0);
      increase+=10;
      
    }

    
    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
    // for(Int_t j = 0; j<2000; j++){
      bch[i]->GetEvent(j);
      charge = 0;
      if(ch[i].selection==0){
        for(Int_t k = 4000/4; k<4400/4.; k++){
          charge += ch[i].wvf[k];
        }
        h[i]->Fill(charge*dtime/sphe);
        hpeak[i]->Fill(ch[i].peak*pow(2,-13));
        h2p->Fill(charge*dtime/sphe,ch[i].peak*pow(2,-13));
        h2->Fill(volts[i],charge*dtime/sphe);
      }
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
    h[i]->Fit(fga[i],"QW");
    fga[i]->Draw("SAME");

    cpeak->cd(i+1);
    hpeak[i]->Draw();
    string histp_name = Form("C2XArapuca A1ch1 Peak - LED %.1f V",volts[i]);
    hpeak[i]->SetNameTitle(histp_name.c_str(),hist_name.c_str());
    hpeak[i]->GetYaxis()->SetTitle("# of events");
    hpeak[i]->GetXaxis()->SetTitle("Amplitude (V)");
    fgap[i] = new TF1(Form("fp%d",i),"gaus(0)",-1,1);
    fgap[i]->SetParameters(0.1,0.1,0.1);
    hpeak[i]->Fit(fgap[i],"QW");
    fgap[i]->Draw("SAME");
  }
  
  vector<Double_t> peak(n);
  vector<Double_t> sigma_peak(n);
  vector<Double_t> Erpeak(n);
  
  vector<Double_t> avgs(n);
  vector<Double_t> sigma_mu(n);
  vector<Double_t> Ersigma_mu(n);
  vector<Double_t> Eravgs(n);
  vector<Double_t> ErSTD(n);
  for(Int_t i = 0; i<n; i++){
    avgs[i] = h[i]->GetMean();
    ErSTD[i] = h[i]->GetStdDev();
    Double_t aux_arror_avgs = h[i]->GetMeanError();
    Double_t aux_arror_std = h[i]->GetStdDevError();

     
    peak[i] = fgap[i]->GetParameter(1);
    sigma_peak[i] = fgap[i]->GetParameter(2);
    Double_t aux_arror_peak = fgap[i]->GetParError(1);
    Double_t aux_arror_std_peak = fgap[i]->GetParError(2);
    if(i==3 || i==4){
      peak[i] = hpeak[i]->GetMean();
      sigma_peak[i] = hpeak[i]->GetStdDev();
      aux_arror_peak = hpeak[i]->GetMeanError();
      aux_arror_std_peak = hpeak[i]->GetStdDevError();

    }
   
    // avgs[i] = fga[i]->GetParameter(1);
    // ErSTD[i] = fga[i]->GetParameter(2);
    // Double_t aux_arror_avgs = fga[i]->GetParError(1);
    // Double_t aux_arror_std = fga[i]->GetParError(2);
    
    
    Eravgs[i] = sqrt(pow(aux_arror_avgs,2)+pow(ErSTD[i],2)/h[i]->GetEntries());
    Erpeak[i] = sqrt(pow(aux_arror_peak,2)+pow(sigma_peak[i],2)/hpeak[i]->GetEntries());

    sigma_mu[i] = 100*ErSTD[i]/avgs[i];
    Ersigma_mu[i] = sigma_mu[i]*sqrt(pow(aux_arror_std/ErSTD[i],2)+pow(aux_arror_avgs/avgs[i],2));
  }
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TCanvas *c2 = new TCanvas("c2","c2",1920,0,1920,1080);

  TGraphErrors *gmeans = new TGraphErrors(n,&volts[0],&avgs[0],0,&Eravgs[0]);
//   TGraphErrors *gSTDs = new TGraphErrors(n,&volts[0],&avgs[0],0,&ErSTD[0]);
  gmeans->GetXaxis()->SetTitle("LED Voltage (V)");
  gmeans->GetYaxis()->SetTitle("Average photo-electrons");
  gmeans->SetMarkerStyle(21);
  gmeans->SetMarkerSize(0.5);
//   gSTDs->SetMarkerSize(0);
//   gSTDs->SetLineWidth(2);
//   gSTDs->Draw("AP []");
  TF1 *fpol1 = new TF1("fpol1","pol1",5,30);
  fpol1->SetParameters(10,10);
  gmeans->Draw("AP");
  gmeans->Fit("fpol1","R");
  
  TCanvas *c3 = new TCanvas("c3","c3",1920,0,1920,1080);

  TGraphErrors *gstdmu = new TGraphErrors(n,&avgs[0],&sigma_mu[0],&Eravgs[0],&Ersigma_mu[0]);
  TF1 *myfunc = new TF1("myfunc","sqrt([0]*[0]+[1]*[1]/x+pow([2]/x,2))",0,200000);
  myfunc->SetParameters(1,2,3);
  gstdmu->GetXaxis()->SetRangeUser(0.001,200000);
  gstdmu->GetYaxis()->SetRangeUser(0,600);
  gstdmu->GetXaxis()->SetTitle("#mu (photo-electrons)");
  gstdmu->GetYaxis()->SetTitle("#sigma/#mu (%)");
  gstdmu->SetMarkerStyle(21);
  gstdmu->SetMarkerSize(0.5);
  gstdmu->Draw("AP");
  gstdmu->Fit("myfunc","R");

  



  TCanvas *c4 = new TCanvas("c4","c4",1920,0,1920,1080);  
  TGraphErrors *gpeakph = new TGraphErrors(n,&avgs[0],&peak[0],&Eravgs[0],&Erpeak[0]);
  TF1 *myfunc2 = new TF1("myfunc2","pol1",0.5,300);
  myfunc2->SetParameters(0.001,0.003);
  // gpeakph->GetXaxis()->SetRangeUser(0.001,200000);
  // gpeakph->GetYaxis()->SetRangeUser(0,600);
  gpeakph->GetXaxis()->SetTitle("#mu (photo-electrons)");
  gpeakph->GetYaxis()->SetTitle("Average amplitude (V)");
  gpeakph->SetMarkerStyle(21);
  gpeakph->SetMarkerSize(0.5);
  gpeakph->Draw("AP");
  gpeakph->Fit("myfunc2","RQ0");
  myfunc2->SetRange(0,250);
  myfunc2->Draw("SAME");

  TCanvas *c2p = new TCanvas("c2p","c2p",1920,0,1920,1080);
  h2p->Draw("colz");
  
  TCanvas *c2c = new TCanvas("c2c","c2c",1920,0,1920,1080);
  h2->Draw("colz");
  
  
}
