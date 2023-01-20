#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


Double_t getstddev(vector<Double_t> vcharge, vector<Double_t> vtime, TF1 *f){
  Double_t result = 0;
  Double_t byhand = 0;
  TH1D *hstd = new TH1D("hstd","hstd",200,0,0);
  Double_t n = 0;
  for(Int_t j = 0; j<vcharge.size(); j++){
    if(vtime[j] >= 200 && vtime[j]<300){
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
void ana_tot(){
  // vector<Int_t> devices = {0,1,2,3};
  vector<Int_t> devices = {0,3};
  // vector<Int_t> devices = {0,3};
  Int_t ndev = devices.size();
  vector<string> files = {"run0_47V00_20ADC_external_trigger_3V3_20ns", "run1_47V00_20ADC_external_trigger_3V5_20ns", "run2_47V00_20ADC_external_trigger_3V7_20ns", "run3_47V00_20ADC_external_trigger_3V9_20ns", "run4_47V00_20ADC_external_trigger_4V1_20ns", "run5_47V00_20ADC_external_trigger_4V3_20ns", "run6_47V00_20ADC_external_trigger_4V5_20ns", "run7_47V00_20ADC_external_trigger_4V7_20ns", "run8_47V00_20ADC_external_trigger_5V_20ns", "run9_47V00_20ADC_external_trigger_6V_20ns", "run10_47V00_20ADC_external_trigger_7V5_20ns", "run11_47V00_20ADC_external_trigger_10V_20ns", "run12_47V00_20ADC_external_trigger_12V5_20ns", "run13_47V00_20ADC_external_trigger_15V_20ns", "run14_47V00_20ADC_external_trigger_17V5_20ns"};

  vector<Double_t> volts = {3.3, 3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 5., 6., 7.5, 10., 12.5, 15.};
  vector<Double_t> saturation_level = {10800,12000,12000,13800};
  vector<Double_t> tot_threshold = {2000,5000,5000,2000};
  // if(mychannel=="Ch5") saturation_level = 12300;
  // if(mychannel=="Ch6") saturation_level = 12600;

  Int_t n = volts.size();

  Double_t dtime = 4;

  string conc = "/analyzed.root";

  vector<ANALYZER*> z(n);
  vector<TH2D*> h2p(ndev);
  ANALYZER *temp = new ANALYZER("temp");

  Double_t sphe = 1; // V.ns
  // vector<Double_t> sphes = {1,1,6064.7,2208.27};
  vector<Double_t> sphes = {1,1,1,1};

  TRandom3 *rd = new TRandom3(1);

  for(Int_t i = 0; i<n; i++){
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z_%.1fV",volts[i]));
    z[i]->setAnalyzer(files[i]);
  }
  for (Int_t i = 0; i < ndev; i++) {
    h2p[i] = new TH2D(Form("h2p_%s",z[0]->schannel[devices[i]].c_str()),Form("h2p_%s",z[0]->schannel[devices[i]].c_str()),2000,-5,4000e3,2000,-100,14000);
  }


  vector<vector<Double_t>> charges(ndev);
  vector<vector<Double_t>> charges_sat(ndev);
  vector<vector<Double_t>> charges_total(ndev);
  vector<vector<Double_t>> tots(ndev);
  vector<vector<Double_t>> tots_sat(ndev);
  vector<vector<Double_t>> tots_total(ndev);
  vector<vector<Double_t>> maxs(ndev);
  vector<vector<Double_t>> maxs_sat(ndev);
  vector<vector<Double_t>> maxs_total(ndev);
  for(Int_t i = 0; i<n; i++){
    Int_t filter = 24;
    Double_t minInt = 10320;
    Double_t maxInt = 10700;
    // if(mychannel == "Ch6"){
    //   minInt = 10320;
    //   maxInt = 10500;
    // }

    Int_t nentries = z[i]->t1->GetEntries();
    nentries = 2000;
    Int_t currentp = charges_total[0].size();
    for(Int_t k = 0; k < ndev; k++){
      Int_t channel = devices[k];
      charges_total[k].resize(charges_total[k].size()+nentries);
      tots_total[k].resize(tots_total[k].size()+nentries);
      maxs_total[k].resize(maxs_total[k].size()+nentries);
      for(Int_t j = 0; j<nentries; j++){
      // if(i == 1 && j==1765) continue;
      // if(i == 1 && j==517) continue;
      // if(i == 12 && j==519) continue;
        sphe = sphes[channel];
        z[i]->getWaveform(j, channel);
        z[i]->applyDenoise(filter);
        z[i]->integrate(channel,minInt,maxInt);
        Double_t tot = z[i]->getTOT(channel,10270,14000,tot_threshold[channel]);
        Double_t charge = z[i]->temp_charge;
        Double_t max = z[i]->temp_max;
        // cout << files[i] << " ev "  << j << " ch " << z[0]->schannel[channel]  << endl;
        // if(i==12 && max < 10000 && z[i]->schannel[channel] == "Ch1") continue;
        h2p[k]->Fill(charge/sphe,max);
        if(max < saturation_level[channel]){
          charges[k].push_back(charge/sphe);
          tots[k].push_back(tot);
          maxs[k].push_back(max);
        }
        else {
          charges_sat[k].push_back(charge/sphe);
          tots_sat[k].push_back(tot);
          maxs_sat[k].push_back(max);
        }
        charges_total[k][currentp+j] = charge/sphe;
        tots_total[k][currentp+j] = tot;
        maxs_total[k][currentp+j] = max;
      }
    }
  }

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  vector<TCanvas *>c2p(ndev);
  vector<TCanvas *> cg(ndev);
  vector<TGraph *> gpeak(ndev);
  vector<TGraph *> gtot(ndev);
  vector<TGraph *> gtot_sat(ndev);
  vector<TF1 *> fexpo(ndev);

  vector<Double_t> corrected_charge = charges_total[0];
  for (Int_t i = 0; i < ndev; i++) {
    string mychannel = z[0]->schannel[devices[i]].c_str();
    c2p[i] = new TCanvas(Form("c2p_%s",mychannel.c_str()),Form("c2p_%s",mychannel.c_str()),1920,0,1920,1080);
    // h2p[i]->GetXaxis()->SetTitle("#mu (photo-electrons)");
    // if(mychannel == "Ch1" || mychannel == "Ch2")
    h2p[i]->GetXaxis()->SetTitle("Charge (ADC*ns)");
    h2p[i]->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
    h2p[i]->Draw("colz");

    cg[i] = new TCanvas(Form("cg_%s",mychannel.c_str()),Form("cg_%s",mychannel.c_str()),1920,0,1920,1080);
    cg[i]->cd();
    // gtot[i] = new TGraph(tots[i].size(),&tots[i][0],&charges[i][0]);
    gtot[i] = new TGraph(tots[i].size(),&tots[i][0],&charges[i][0]);
    gtot_sat[i] = new TGraph(tots_sat[i].size(),&tots_sat[i][0],&charges_sat[i][0]);

    gtot[i]->SetMarkerStyle(21);
    gtot[i]->SetMarkerSize(0.7);
    gtot_sat[i]->SetMarkerStyle(21);
    gtot_sat[i]->SetMarkerSize(0.7);
    gtot_sat[i]->SetMarkerColor(kRed);
    // gtot[i] = new TGraph(tots[i].size(),&charges[i][0],&maxs[i][0]);
    gtot[i]->Draw("AP");
    gtot_sat[i]->Draw("SAME P");

    fexpo[i] = new TF1(Form("fexpo_%d",i),"expo",100,340);
    gtot[i]->Fit(Form("fexpo_%d",i),"RQ");
    if (i == 0) {
      fexpo[0]->SetRange(50,1600);
      fexpo[0]->SetNpx(400);
      Double_t stddev = getstddev(charges[0],tots[0],fexpo[0]);
      // I am assuming the std is not fixed, but proportional to the amplitude of the charge. So:
      cout << fexpo[0]->Eval(200) << endl;
      cout << stddev << endl;
      stddev = stddev/fexpo[0]->Eval(200);
      Double_t refval = 0;
      for(Int_t j = 0; j<charges_total[0].size(); j++){
        if(maxs_total[0][j] > saturation_level[0])
        {
          if (tots_total[0][j]!=0) {
            refval = fexpo[0]->Eval(tots_total[0][j]);
            corrected_charge[j] = rd->Gaus(refval,stddev*refval);
          }
        }
      }
      TGraph *gtot_corr = new TGraph(charges_total[0].size(), &tots_total[0][0], &corrected_charge[0]);
      gtot_corr->SetMarkerStyle(21);
      gtot_corr->SetMarkerSize(0.7);
      gtot_corr->SetMarkerColor(kBlue);
      gtot_corr->Draw("SAME P");

    }


  }

  if(ndev>1){
    TMultiGraph *gm = new TMultiGraph("gm","gm");


    TCanvas *ccomp = new TCanvas("ccomp","ccomp");
    TGraph *gcomp = new TGraph(charges_total[0].size(),&charges_total[1][0], &charges_total[0][0]);
    gcomp->GetXaxis()->SetTitle("Charge Ch6 (ADC*ns)");
    gcomp->GetYaxis()->SetTitle("Charge Ch1 (ADC*ns)");
    gm->Add(gcomp);
    // gcomp->Draw("AP");


    TGraph *gcomp_corr = new TGraph(charges_total[0].size(), &charges_total[1][0], &corrected_charge[0]);
    gcomp_corr->SetMarkerColor(kBlue);
    gm->Add(gcomp_corr);
    // gcomp_corr->Draw("SAME P");

    gm->Draw("AP");


  }


}
