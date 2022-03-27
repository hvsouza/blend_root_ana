#define memorydepth 1252
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


class OUTPUT{
public:
  Double_t charge;
  Double_t slow;
    
  OUTPUT();
  const char *tobranch = Form("charge/D:slow/D");
};

OUTPUT::OUTPUT(){
  charge = 0;
  slow = 0;
  }


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

Double_t getstddev(vector<Double_t> vcharge, vector<Double_t> vtime, TF1 *f){
  Double_t result = 0;
  Double_t byhand = 0;
  TH1D *hstd = new TH1D("hstd","hstd",50,-100,100);
  Double_t n = 0;
  for(Int_t j = 0; j<vcharge.size(); j++){
    if(vtime[j] >= 70 && vtime[j]<250){
      hstd->Fill(vcharge[j]-(f->Eval(vtime[j])));
      // cout << vcharge[j] << " " << f->Eval(vtime[j]) << " " << vtime[j] << endl;
      byhand+=(vcharge[j]-(f->Eval(vtime[j])))*(vcharge[j]-(f->Eval(vtime[j])));
      n+=1;
    }
  }
  byhand = sqrt(byhand/n);
  // TCanvas *ctemp = new TCanvas("ctemp");  
  // hstd->Draw();
  // ctemp->Print(Form("testing_%2f.root",vcharge[0]));
  result = hstd->GetStdDev();
  delete hstd;
  // delete ctemp;


  // cout << result << " " << byhand << endl;
  return result;
 
  

}


void evaluate_undershoot_amplitude(){
  // vector<Double_t> volts = {10,12.5,15,17.5,20,22,25,27.5,30,5,7.5};
  vector<Double_t> volts = {0,10};
  Int_t n = volts.size();
  TFile *fout = new TFile("C2_output_charges.root","RECREATE");
  TTree *tout = new TTree("t1","t1");
  vector<OUTPUT> out(n);
  vector<TBranch*> br(n);
  for (Int_t i = 0; i<n; i++) br[i] = tout->Branch(Form("E%.0f",volts[i]),&out[i],out[0].tobranch);
  // vector<string> filesC2 = {"C1XArapuca_Efield_0kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_10kV_A4ch2_250MHz_cosmic00000","C1XArapuca_Efield_5kV_A4ch2_250MHz_cosmic00000"};
  vector<string> filesC2 = {"C2XArapuca_Efield_0kV_A1ch1_250MHz_cosmic00000","C2XArapuca_Efield_10kV_A1ch1_250MHz_cosmic00000","C2XArapuca_Efield_5kV_A1ch1_250MHz_cosmic00000"};
  string conc = "/analyzed.root";

  vector<TFile*> f(n);
  vector<TTree*> t1(n);
  vector<ADC_DATA> ch(n);
  vector<TBranch *> bch(n);
  vector<TH1D*> h(n);
  vector<TH1D*> hcor(n);
  vector<TGraph*> g1(n);
  vector<TGraph*> g1_sat(n);
  vector<TGraph*> g1_corrected(n);
  vector<TGraph*> g2(n);
  vector<TGraph*> g3(n);
  vector<TGraph*> g3_sat(n);
  vector<TMultiGraph*> gunders(n);
  vector<TMultiGraph*> grise(n);
  vector<TF1*> fpol1(n);
  vector<TF1*> fexpo(n);
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
  vector<vector<Double_t>> vcharge_normal(n);
  vector<vector<Double_t>> vcharge_sat(n);
  vector<vector<Double_t>> vcharge_corrected(n);
  vector<vector<Double_t>> vundercharge_normal(n);
  vector<vector<Double_t>> vundercharge_sat(n);
  vector<vector<Double_t>> vpeak_normal(n);
  vector<vector<Double_t>> vpeak_sat(n);
  vector<vector<Double_t>> vunder_normal(n);
  vector<vector<Double_t>> vunder_sat(n);
  vector<vector<Double_t>> time_rise_normal(n);
  vector<vector<Double_t>> time_rise_sat(n);
  Double_t time_trigger = 0.4;
  Double_t undershoot = 0;

  TRandom3 *rd = new TRandom3(1);


  for(Int_t i = 0; i<n; i++){
    filesC2[i] = filesC2[i]+conc;
    f[i] = new TFile(filesC2[i].c_str(),"READ");
    t1[i] = (TTree*)f[i]->Get("t1");
    bch[i] = t1[i]->GetBranch("Ch1");
    bch[i]->SetAddress(&ch[i]);
    h[i] = new TH1D(Form("h%d",i),Form("h%d",i),50,0,800);
    hcor[i] = new TH1D(Form("hcor%d",i),Form("hcor%d",i),100,0,4000);
    // h[i] = new TH1D(Form("h%d",i),Form("h%d",i),20000,-1,1);
    tout->SetEntries(t1[i]->GetEntries());
    for(Int_t j = 0; j<t1[i]->GetEntries(); j++){

    // for(Int_t j = 1; j<3; j++){
      bch[i]->GetEvent(j);
      charge = 0;
      charge_undershoot = 0;
      undershoot = 1e12;
      Int_t never_negative = true;
      Double_t newBase = getbaseline(ch[i].wvf);
      for(Int_t k = 800/4; k<2000/4.; k++){
      // for(Int_t k = 840/4; k<3600/4.; k++){ \\ Ajib

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

      Bool_t triggered = false;
      Double_t time_mark = 0;
      for(Int_t k = 700/4; k<2400/4; k++){
	if(ch[i].wvf[k]>time_trigger && triggered == false){
	  triggered = true;
          time_mark = k*4;
	  // break;
	}
        else if(triggered==true && ch[i].wvf[k]<time_trigger && (k*4-time_mark)>=12){
          time_mark = k*4 - time_mark;

          break;
        }
      }
      if(ch[i].peak<=1.4){
        if(triggered==false) time_rise_normal[i].push_back(0);
        else time_rise_normal[i].push_back(time_mark);
        vcharge_normal[i].push_back(charge*4);
        vpeak_normal[i].push_back(ch[i].peak);
        vunder_normal[i].push_back(undershoot);
        vundercharge_normal[i].push_back(charge_undershoot*4);
        h[i]->Fill(charge_undershoot*4);
        hcor[i]->Fill(charge*4);
        // filling normally here
        out[i].charge = charge*4;
        out[i].slow = charge_undershoot*4;
        br[i]->Fill();
          
      }
      else{
        if(triggered==false) time_rise_sat[i].push_back(0);
        else time_rise_sat[i].push_back(time_mark);
        vcharge_sat[i].push_back(charge*4);
        vpeak_sat[i].push_back(ch[i].peak);
        vunder_sat[i].push_back(undershoot);
        vundercharge_sat[i].push_back(charge_undershoot*4);
        h[i]->Fill(charge_undershoot*4);
      }
    }
    



    // c1->cd(i+1);
    c[aux+i] = new TCanvas();
    c[aux+i]->cd();
    aux++;
    gunders[i] = new TMultiGraph(Form("Undershoot %s",filesC2[i].c_str()),Form("Undershoot %s",filesC2[i].c_str()));
    // g3[i] = new TGraph(vcharge[i].size(),&vundercharge[i][0], &vcharge[i][0]);
    g3[i] = new TGraph(vcharge_normal[i].size(), &vcharge_normal[i][0],&vundercharge_normal[i][0]);
    g3_sat[i] = new TGraph(vcharge_sat[i].size(), &vcharge_sat[i][0],&vundercharge_sat[i][0]);
    fpol1[i] = new TF1(Form("f_%.0f_kV",volts[i]),"pol1",40,1000);
    
    g3[i]->SetNameTitle(Form("%.0f kV",volts[i]),Form("%.0f kV",volts[i]));
    g3[i]->SetMarkerStyle(21);
    g3[i]->SetMarkerSize(0.7);
    g3_sat[i]->SetMarkerStyle(21);
    g3_sat[i]->SetMarkerSize(0.7);
    g3_sat[i]->SetMarkerColor(kRed);
    // g3[i]->Draw("AP");
    // g3_sat[i]->Draw("SAME P");
    g3[i]->Fit(Form("f_%.0f_kV",volts[i]),"RQ");
    gunders[i]->Add(g3[i]);
    gunders[i]->Add(g3_sat[i]);

    gunders[i]->GetXaxis()->SetTitle("Charge (V.ns)");
    gunders[i]->GetYaxis()->SetTitle("Slow component charge (V.ns)");
    gunders[i]->Draw("AP");

    
    c[aux+i] = new TCanvas();
    c[aux+i]->cd();
    grise[i] = new TMultiGraph(Form("TOT %s",filesC2[i].c_str()),Form("TOT %s",filesC2[i].c_str()));
    g1[i] = new TGraph(vcharge_normal[i].size(),&time_rise_normal[i][0], &vcharge_normal[i][0]);
    g1_sat[i] = new TGraph(vcharge_sat[i].size(),&time_rise_sat[i][0], &vcharge_sat[i][0]);
    g1[i]->SetMarkerStyle(21);
    g1[i]->SetMarkerSize(0.7);
    g1_sat[i]->SetMarkerStyle(21);
    g1_sat[i]->SetMarkerSize(0.7);
    g1_sat[i]->SetMarkerColor(kRed);
    // g1[i]->Draw("AP");
    // g1_sat[i]->Draw("SAME P");
    grise[i]->Add(g1[i]);
    grise[i]->Add(g1_sat[i]);
    fexpo[i] = new TF1(Form("fexpo_%.0fkV",volts[i]),"expo",70,250);
    g1[i]->Fit(Form("fexpo_%.0fkV",volts[i]),"RQ");
    Double_t stddev = getstddev(vcharge_normal[i],time_rise_normal[i],fexpo[i]);
    // I am assuming the std is not fixed, but proportional to the amplitude of the charge. So:
    stddev = stddev/fexpo[i]->Eval(200);
    Double_t refval = 0;
    for(Int_t j = 0; j<vcharge_sat[i].size(); j++){
      refval = fexpo[i]->Eval(time_rise_sat[i][j]); 
      vcharge_corrected[i].push_back(rd->Gaus(refval,stddev*refval));
      hcor[i]->Fill(vcharge_corrected[i][j]);

      out[i].charge = vcharge_corrected[i][j];
      out[i].slow = vundercharge_sat[i][j];
      br[i]->Fill();
    }
        
    g1_corrected[i] = new TGraph(vcharge_sat[i].size(),&time_rise_sat[i][0], &vcharge_corrected[i][0]);
    g1_corrected[i]->SetMarkerStyle(21);
    g1_corrected[i]->SetMarkerSize(0.7);
    g1_corrected[i]->SetMarkerColor(kBlue);
    grise[i]->Add(g1_corrected[i]);
    gPad->SetLogy(1);
    grise[i]->GetYaxis()->SetRangeUser(4,5000);
    grise[i]->GetXaxis()->SetTitle("Time over threshold (ns)");
    grise[i]->GetYaxis()->SetTitle("Charge (V.ns)");
    grise[i]->Draw("AP");
    aux++;


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
  h[0]->GetXaxis()->SetTitle("Slow compoment charge (V*ns)");
  h[0]->GetYaxis()->SetTitle("# of events");
  h[0]->SetLineColor(kRed);
  h[1]->SetLineColor(kBlack);
  h[0]->SetLineWidth(2);
  h[1]->SetLineWidth(2);
  // h[2]->Draw("same");


  TCanvas *c8 = new TCanvas();
  hcor[0]->Draw();
  hcor[1]->Draw("same");
  hcor[0]->GetXaxis()->SetTitle("Charge (V*ns)");
  hcor[0]->GetYaxis()->SetTitle("# of events");
  hcor[0]->SetLineColor(kRed);
  hcor[1]->SetLineColor(kBlack);
  hcor[0]->SetLineWidth(2);
  hcor[1]->SetLineWidth(2);  
  

  fout->WriteObject(tout,"t1","TObject::KOverwrite");
  
}
