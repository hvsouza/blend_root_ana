#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"



Double_t getbaseline(Double_t v[]){
  Double_t result = 0;
  TH1D *hbase = new TH1D("hbase","hbase",2000,-0.1,0.1);
  for(Int_t i=0; i<800/4; i++) hbase->Fill(v[i]);
  Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());  
  Int_t bins=0;
  for(Int_t i=0; i<800/4;){
    if(v[i] > res0 + 0.02) {
      i+=100/4;
    }
    else{
      result += v[i];
      bins++;
      i++;
    }
  }
  result/=bins;
  if(bins > (800/4)/2){
    delete hbase;
    return result;
  }
  else{
    delete hbase;    
    return res0;
  }

}

void evaluate_undershoot_amplitude(){
  // vector<Double_t> volts = {10,12.5,15,17.5,20,22,25,27.5,30,5,7.5};
  vector<Double_t> volts = {0,10};
  Int_t n = volts.size();
  vector<string> filesC2 = {"C1XArapuca_Efield_0kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_10kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_5kV_A4ch2_250MHz_cosmic00000"};
//   vector<string> filesC2 = {"C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl10V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl12p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl15V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl17p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl20V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl22p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl25V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl27p5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl30V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl5V00000","C1XArapuca_Efield_0kV_CRPon_A1ch1_250MHz_LEDwith20nsampl7p5V00000"};
  // vector<string> filesC2 = {"C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl10V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl12p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl15V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl17p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl20V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl22p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl25V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl27p5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl30V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl5V00000","C2XArapuca_Efield_0kV_CRPon_A4ch2_250MHz_LEDwith20nsampl7p5V00000"};

  string conc = "/analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);
  vector<TH1D*> h(n);
  vector<TGraph*> g1(n);
  vector<TGraph*> g2(n);
  vector<TGraph*> g3(n);
  vector<TF1*> fpol1(n);
  // gStyle->SetPadTickX(0);
  // gStyle->SetPadTickY(0);
  // TCanvas *c1 = new TCanvas();
  vector<TCanvas*> c(n*2);
  vector<TF1*> fga(n);
  Double_t charge = 0;
  Double_t charge_undershoot = 0;
  Int_t m1=0,m2=0;
  m1=m2=n/2;
  // cout << m2 << endl;
  if(n%2==1) m1+=1;
  if (m1*m2 < n) m2+=1;
  
  // c1->Divide(m1,m2);
  // c1->Update();
  Int_t aux = 0;
  vector<vector<Double_t>> vcharge(n);
  vector<vector<Double_t>> vundercharge(n);
  vector<vector<Double_t>> vpeak(n);
  vector<vector<Double_t>> vunder(n);
  vector<vector<Double_t>> time_rise(n);
  Double_t time_trigger = 0.1;
  Double_t undershoot = 0;
  for(Int_t i = 0; i<n; i++){
    filesC2[i] = filesC2[i]+conc;
    f[i] = new TFile(filesC2[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch("Ch1");
    bch[i]->SetAddress(&ch[i]);
    h[i] = new TH1D(Form("h%d",i),Form("h%d",i),200,0,800);
    // h[i] = new TH1D(Form("h%d",i),Form("h%d",i),20000,-1,1);
    
    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){
      bch[i]->GetEvent(j);
      charge = 0;
      charge_undershoot = 0;
      undershoot = 1e12;
      Int_t never_negative = true;
      Double_t newBase = getbaseline(ch[i].wvf);
      for(Int_t k = 900/4; k<2000/4.; k++){

	charge += ch[i].wvf[k];
	if(ch[i].wvf[k]<=undershoot && k>1100/4) undershoot = ch[i].wvf[k];
	if(k>1150/4){
	  // if(i==0)
	    charge_undershoot+=ch[i].wvf[k];
	  // else
	    // charge_undershoot+=ch[i].wvf[k];
	}
     
	// if(ch[i].wvf[k]<newBase && k>1050/4){
	  // never_negative = false;
	// }
        // if(never_negative==true) charge += ch[i].wvf[k];
	// else charge_undershoot+=ch[i].wvf[k];
      }
      // if(ch[i].peak>1.3) continue;
      Bool_t triggered = false;
      for(Int_t k = 700/4; k<1200; k++){
	if(ch[i].wvf[k]>time_trigger){
	  time_rise[i].push_back(k*4);
	  triggered = true;
	  break;
	}
      }
      if(triggered==false) time_rise[i].push_back(0);
      vcharge[i].push_back(charge*4);
      vpeak[i].push_back(ch[i].peak);
      vunder[i].push_back(undershoot);
      vundercharge[i].push_back(charge_undershoot*4);
      h[i]->Fill(charge_undershoot*4);
    }



    // c1->cd(i+1);
    c[aux+i] = new TCanvas();
    c[aux+i]->cd();
    aux++;
    // g3[i] = new TGraph(vcharge[i].size(),&vundercharge[i][0], &vcharge[i][0]);
    g3[i] = new TGraph(vcharge[i].size(), &vcharge[i][0],&vundercharge[i][0]);
    fpol1[i] = new TF1(Form("f_%.0f_kV",volts[i]),"pol1",40,1000);
    
    g3[i]->SetNameTitle(Form("%.0f kV",volts[i]),Form("%.0f kV",volts[i]));
    g3[i]->GetXaxis()->SetTitle("Charge (V.ns)");
    g3[i]->GetYaxis()->SetTitle("Slow component charge (V.ns)");
    g3[i]->Draw("AP*");
    g3[i]->Fit(Form("f_%.0f_kV",volts[i]),"RQ");
    
    // c[aux+i] = new TCanvas();
    // c[aux+i]->cd();
    // g1[i] = new TGraph(vcharge[i].size(),&time_rise[i][0], &vcharge[i][0]);
    // g1[i]->GetXaxis()->SetTitle("Rising edge (ns)");
    // g1[i]->GetYaxis()->SetTitle("Charge (V.ns)");
    // g1[i]->Draw("AP*");
    // aux++;


    // c[aux+i] = new TCanvas();
    // c[aux+i]->cd();
    // g2[i] = new TGraph(vcharge[i].size(),&vcharge[i][0], &vpeak[i][0]);
    // g2[i]->GetXaxis()->SetTitle("Charge (V.ns)");
    // g2[i]->GetYaxis()->SetTitle("Max. amplitude (V)");
    // g2[i]->Draw("AP*");
    // aux++;
   
    
    // c[aux+i] = new TCanvas();
    // g2[i] = new TGraph(vcharge[i].size(),&vunder[i][0], &vcharge[i][0]);
    // g2[i]->GetYaxis()->SetTitle("Charge (V.ns)");
    // g2[i]->GetXaxis()->SetTitle("Undershoot amplitude (V)");
    // g2[i]->Draw("AP*");

    // aux++;
   
//     h[i]->GetXaxis()->SetRangeUser(-0.04,0.04);
    // h[i]->SetNameTitle(Form("X-ARAPUCA C1 - LED %1.f V",volts[i]),Form("X-ARAPUCA C1 - LED %1.f V",volts[i]));
    // h[i]->GetYaxis()->SetTitle("# of events");
    // h[i]->GetXaxis()->SetTitle("Amplitude (V)");
    // fga[i] = new TF1(Form("f%d",i),"gaus(0)",-1,1);
    // fga[i]->SetParameters(0.1,0.1,0.1);
    // h[i]->Fit(fga[i]);
    // fga[i]->Draw("SAME");
  }

  TCanvas *c7 = new TCanvas();
  h[0]->Draw();
  h[1]->Draw("same");
  // h[2]->Draw("same");

  
  
  
  
}
