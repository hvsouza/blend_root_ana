#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void create_histograms_argon2x2(){
  gStyle->SetCanvasPreferGL(kFALSE);
  bool just_a_test = false;
  string mychannel = "Ch0";
  Double_t minInt = 10380;
  Double_t maxInt = 11400;
  if(mychannel == "Ch0"){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch2"){minInt = 10280; maxInt = 10700;}
  if(mychannel == "Ch3"){minInt = 10320; maxInt = 10500;}
  if(mychannel == "Ch4"){minInt = 10300; maxInt = 10720;}
  vector<string> devices = {
                            "miniArapuca 37V (A1ch1)",
                            "miniArapuca 37V (Argon4)",
                            "xArapuca V1 36V (DCemVD1.0)",
                            "xArapuca V4 Ch2 (DCemArgon4+INFN)",
                            "xArapuca V4 Ch1 (DCemArgon4+INFN)",
                            "xArapuca V5 Ch2 (DCem1.2+LBL)",
                            "xArapuca V5 Ch1 (DCem1.2+LBL)",
                            "flex 46V (DCemVD1.2)",};

  vector<Double_t> saturations = {12000, 12000, 12600, 12600, 12000, 12000, 12000, 100};
  vector<Double_t> sphes = {7100,400, 7100, 1, 5878, 0, 0, 0};

  vector<string> files = {"run23_all_devices_275nm_20ns_4V0", "run24_all_devices_275nm_20ns_4V5", "run25_all_devices_275nm_20ns_5V0", "run26_all_devices_275nm_20ns_5V5", "run27_all_devices_275nm_20ns_6V0", "run28_all_devices_275nm_20ns_6V5", "run29_all_devices_275nm_20ns_7V0", "run30_all_devices_275nm_20ns_7V5", "run31_all_devices_275nm_20ns_8V0", "run32_all_devices_275nm_20ns_8V5", "run33_all_devices_275nm_20ns_9V0", "run34_all_devices_275nm_20ns_10V0", "run35_all_devices_275nm_20ns_12V0", "run36_all_devices_275nm_20ns_14V0", "run37_all_devices_275nm_20ns_16V0", "run38_all_devices_275nm_20ns_18V0", "run39_all_devices_275nm_20ns_20V0", "run40_all_devices_275nm_20ns_22V5", "run41_all_devices_275nm_20ns_25V0", "run42_all_devices_275nm_20ns_27V5", "run43_all_devices_275nm_20ns_30V0",};
  vector<Double_t> volts = {4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.5, 25.0, 27.5, 30.0};
  Double_t saturation_level;

  Int_t n = volts.size();

  // vector<string> files(n);
  // for(Int_t i=0; i<n; i++){
  //   if(static_cast<Int_t>(volts[i])-volts[i]!=0) files[i] = Form("%.1fV",volts[i]);
  //   else files[i] = Form("%.0fV",volts[i]);
  // }
  vector<Double_t> ranges = {-4000,4000};

  Int_t nbins;
  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);

  vector<TH1D*> h(n);
  vector<TH1D*> hpeak(n);
  TH2D *h2p = new TH2D("h2p","h2p",4000,-5,4000,2000,-100,14000);
  TH2D *h2 = new TH2D("h2","h2",300,0,30,8000,-50,4000);
  TGraphErrors *gpeaks = new TGraphErrors();
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  TCanvas *cpeak = new TCanvas("cpeak","cpeak",1920,0,1920,1080);
  vector<TF1*> fga(n);
  vector<TF1*> fgap(n);
  Double_t charge = 0;
  Double_t max = 0;
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


  Double_t sphe = 6046; // V.ns
  Double_t increase = 0;
  Int_t nbins_d = 100;
  Int_t gpts = 0;
  ANALYZER a;
  for(Int_t i = 0; i<n; i++){
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z%.2f",volts[i]));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    Int_t kch = z[i]->kch;
    saturation_level = saturations[z[i]->getIdx()];
    sphe = sphes[z[i]->getIdx()];

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

    DENOISE dn;
    Int_t filter = 16;
    vector<Double_t> wvf(memorydepth);
    Int_t nentries = (just_a_test) ? 1000 : z[i]->nentries;
    for(Int_t j = 0; j<nentries; j++){
      z[i]->getWaveform(j,kch);
      charge = 0;
      if(z[i]->ch[kch].selection==0){

        z[i]->applyDenoise(filter);
        z[i]->integrate(minInt,maxInt);
        charge = z[i]->temp_charge;
        max = z[i]->temp_max;
        h[i]->Fill(charge/sphe);
        hpeak[i]->Fill(max);
        h2p->Fill(charge/sphe,max);
        gpeaks->SetPoint(gpts,charge/sphe,max);
        gpeaks->SetPointError(gpts,sqrt(charge/sphe),sqrt(max));
        // cout << gpts << " " << gpeaks->GetN() << " " << h2p->GetEntries() << endl;
        gpts++;
        h2->Fill(volts[i],charge/sphe);
      }
    }
    c1->cd(i+1);
    h[i]->Draw();
//     h[i]->GetXaxis()->SetRangeUser(-0.04,0.04);
    string hist_name = Form("%s - LED %.1f V",devices[z[i]->getIdx()].c_str(),volts[i]);
    h[i]->SetNameTitle(hist_name.c_str(),hist_name.c_str());
    h[i]->GetYaxis()->SetTitle("# of events");
    h[i]->GetXaxis()->SetTitle("Photo-electrons");
    fga[i] = new TF1(Form("f%d",i),"gaus(0)",-1,1);
    fga[i]->SetParameters(0.1,0.1,0.1);
    // h[i]->Fit(fga[i],"QW");
    // fga[i]->Draw("SAME");

    cpeak->cd(i+1);
    hpeak[i]->Draw();
    string histp_name = Form("%s peak - LED %.1f V",devices[z[i]->getIdx()].c_str(),volts[i]);
    hpeak[i]->SetNameTitle(histp_name.c_str(),hist_name.c_str());
    hpeak[i]->GetYaxis()->SetTitle("# of events");
    hpeak[i]->GetXaxis()->SetTitle("Amplitude (ADC Channels)");
    fgap[i] = new TF1(Form("fp%d",i),"gaus(0)",-1,1);
    fgap[i]->SetParameters(0.1,0.1,0.1);
    // hpeak[i]->Fit(fgap[i],"QW");
    // fgap[i]->Draw("SAME");
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


    // peak[i] = fgap[i]->GetParameter(1);
    // sigma_peak[i] = fgap[i]->GetParameter(2);
    // Double_t aux_arror_peak = fgap[i]->GetParError(1);
    // Double_t aux_arror_std_peak = fgap[i]->GetParError(2);
    // if(i==3 || i==4){
      peak[i] = hpeak[i]->GetMean();
      sigma_peak[i] = hpeak[i]->GetStdDev();
      Double_t aux_arror_peak = hpeak[i]->GetMeanError();
      Double_t aux_arror_std_peak = hpeak[i]->GetStdDevError();

    // }

    // avgs[i] = fga[i]->GetParameter(1);
    // ErSTD[i] = fga[i]->GetParameter(2);
    // Double_t aux_arror_avgs = fga[i]->GetParError(1);
    // Double_t aux_arror_std = fga[i]->GetParError(2);


    // Eravgs[i] = sqrt(pow(aux_arror_avgs,2)+pow(ErSTD[i],2)/* /h[i]->GetEntries() */);
    Eravgs[i] = ErSTD[i];
    // Eravgs[i] = ErSTD[i]/sqrt(avgs[i]);
    // Erpeak[i] = sqrt(pow(aux_arror_peak,2)+pow(sigma_peak[i],2)/hpeak[i]->GetEntries());
    Erpeak[i] = sigma_peak[i];

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
  myfunc->SetParLimits(0,0,1e3);
  myfunc->SetParLimits(1,0,1e3);
  myfunc->SetParLimits(2,0,1e3);
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
  gpeakph->GetYaxis()->SetTitle("Average amplitude (ADC Channels)");
  gpeakph->SetMarkerStyle(21);
  gpeakph->SetMarkerSize(0.5);
  gpeakph->Draw("AP");
  gpeakph->GetXaxis()->SetLimits(0.01,3000);
  // gpeakph->Fit("myfunc2","RQ0");
  // myfunc2->SetRange(0,250);
  // myfunc2->Draw("SAME");

  TCanvas *c2p = new TCanvas("c2p","c2p",1920,0,1920,1080);
  h2p->GetXaxis()->SetTitle("#mu (photo-electrons)");
  h2p->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
  h2p->Draw("colz");

  TCanvas *c2c = new TCanvas("c2c","c2c",1920,0,1920,1080);
  h2->GetXaxis()->SetTitle("LED voltage (V)");
  h2->GetYaxis()->SetTitle("Photo-electrons");
  h2->Draw("colz");

  c2p->cd();
  gpeaks->Draw("L SAME");
  gpeaks->SetLineWidth(0);
  TF1 *myfunc3 = new TF1("myfunc3","pol1",-20,300);
  gpeaks->Fit("myfunc3");


}
