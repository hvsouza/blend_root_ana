#define memorydepth 5000
// #define DUNESTYLE_ENABLE_AUTOMATICALLY 0
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

bool just_a_test = false;
Double_t threshold = 300;
Double_t minFit = 100; // changed later
Double_t rebin = 8;
// Double_t rebin = 8; // lower than 600
Double_t filter = 16;
// gStyle->SetCanvasPreferGL(kFALSE);
string mychannel = "Ch3";
vector<string> devices = {"0","0","miniArapuca 37V (A1ch1)","miniArapuca 47V (Argon4)", "xArapuca v4 ch1 (DCemArgon4)", "xArapuca v4 ch2 (DCemArgon4)", "xArapuca v5 ch1 (DCemSimp3)", "xArapuca v5 ch2 (DCemSimp3)"};
vector<Double_t> saturations = {0, 0, 12600, 2000, 12000, 12000, 12000, 100};
vector<Double_t> sphes = {0,0, 7100, 3.8, 441*6, 0, 0,};

vector<string> files = {"run0_all_devices_led_365nm_20ns_3V15", "run1_all_devices_led_365nm_20ns_3V20", "run2_all_devices_led_365nm_20ns_3V30", "run3_all_devices_led_365nm_20ns_3V50", "run4_all_devices_led_365nm_20ns_3V70", "run5_all_devices_led_365nm_20ns_3V90", "run6_all_devices_led_365nm_20ns_4V10", "run7_all_devices_led_365nm_20ns_4V30", "run8_all_devices_led_365nm_20ns_4V50", "run9_all_devices_led_365nm_20ns_4V70", "run10_all_devices_led_365nm_20ns_4V90", "run11_all_devices_led_365nm_20ns_5V20", "run12_all_devices_led_365nm_20ns_6V20", "run13_all_devices_led_365nm_20ns_7V50", "run14_all_devices_led_365nm_20ns_10V00", "run15_all_devices_led_365nm_20ns_12V50", "run16_all_devices_led_365nm_20ns_15V00", "run17_all_devices_led_365nm_20ns_17V50", "run18_all_devices_led_365nm_20ns_23V00", "run19_all_devices_led_365nm_20ns_30V00"};
vector<Double_t> volts = {3.15, 3.20, 3.30, 3.50, 3.70, 3.90, 4.10, 4.30, 4.50, 4.70, 4.90, 5.20, 6.20, 7.50, 10.00, 12.50, 15.00, 17.50, 23.00, 30.00};
Double_t saturation_level;

Int_t n = volts.size();

string conc = "/analyzed.root";

vector<ANALYZER*> z(n);

Int_t bins_pe = 1000;
Int_t bins_tot = 1000;
Double_t minpe = 0;
Double_t mintot = 0;

Double_t maxpe = 3000;
Double_t maxtot = 4000;
TH2D *htot = new TH2D("htot","htot", bins_pe, minpe, maxpe, bins_tot, mintot, maxtot);
TH2D *htot_sat = new TH2D("htot_sat","htot_sat", bins_pe, minpe, maxpe, bins_tot, mintot, maxtot);
TH2D *htot_before_sat = new TH2D("htot_before_sat","htot_before_sat", bins_pe, minpe, maxpe, bins_tot, mintot, maxtot);


TH2D *hQ = new TH2D("hQ","hQ",bins_tot, mintot, maxtot, bins_pe, minpe, maxpe);
TH2D *hQ_sat = new TH2D("hQ_sat","hQ_sat",bins_tot, mintot, maxtot, bins_pe, minpe, maxpe);
TH2D *hQ_before_sat = new TH2D("hQ_before_sat","hQ_before_sat",bins_tot, mintot, maxtot, bins_pe, minpe, maxpe);

TH2D *hshow_amp_vs_pe = new TH2D("hshow_amp_vs_pe","hshow_amp_vs_pe",4000,-5,4000,2000,-100,14000);

Int_t nx = htot->GetNbinsX();
Int_t ny = htot->GetNbinsY();
Double_t hmin = htot->GetXaxis()->GetBinLowEdge(1);
Double_t hmax = htot->GetXaxis()->GetBinLowEdge(nx) + htot->GetXaxis()->GetBinWidth(nx);
TRandom3 *r3 = new TRandom3(3);

TH1D *htot1 = new TH1D("htot1","htot1",nx/rebin,hmin,hmax);
TH1D *htot1_sat = new TH1D("htot1_sat","htot1_sat",nx/rebin,hmin,hmax);
TH1D *htot1_before_sat = new TH1D("htot1_before_sat","htot1_before_sat",nx/rebin,hmin,hmax);

vector<vector<Double_t>> h1tot(2);
vector<vector<Double_t>> h1tot_sat(2);
vector<vector<Double_t>> h1tot_before_sat(2);
vector<vector<Double_t>> hpe(2);
vector<vector<Double_t>> hpe_sat(2);
vector<vector<Double_t>> hpe_before_sat(2);
TH1D *hdev = new TH1D("hdev","hdev",200,-1,1);
TH1D *hdev_sat = new TH1D("hdev_sat","hdev_sat",200,-1,1);
ANALYZER a("a");


Double_t minInt = 0;
Double_t maxInt = 0;

vector<Double_t> gpe;
vector<Double_t> er_gpe;
vector<Double_t> gtot;
vector<Double_t> er_gtot;

vector<Double_t> gpe_sat;
vector<Double_t> gpe_sat_r;
vector<Double_t> er_gpe_sat;
vector<Double_t> gtot_sat;
vector<Double_t> er_gtot_sat;

vector<Double_t> gpe_corr;
vector<Double_t> gpe_corr_r;
vector<Double_t> er_gpe_corr;
vector<Double_t> gtot_corr;
vector<Double_t> er_gtot_corr;
Double_t sphe = 6046; // V.ns
void setup_histograms(TH2D *htot, TH1D *htot1, vector<vector<Double_t>> &_h1tot, vector<vector<Double_t>> &_hpe){
  Int_t nx = htot->GetNbinsX();
  Int_t ny = htot->GetNbinsY();
  Double_t sum = 0;
  Double_t div = 0;
  Double_t stddev = 0;
  Double_t n = 0;
  Int_t nrebins = 0;
  int ndiv = 0;
  for(Int_t i = 0; i < nx; i++) { // for loop in the X axis of the histogram;
    for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
      Double_t weight = htot->GetBinContent(i+1,j+1);
      if (weight != 0){
        Double_t bnpes = htot->GetYaxis()->GetBinCenter(j);
        sum += weight*bnpes;
        div += weight;
      }
    }
    ndiv++;
    if(ndiv == rebin){
      if(div == 0){sum = 0; div = 1;}
      for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
        Double_t weight = htot->GetBinContent(i+1,j+1);
        if (weight != 0){
          Double_t bnpes = htot->GetYaxis()->GetBinCenter(j);
          stddev += weight*(bnpes-sum/div)*(bnpes-sum/div);
          n += 1;
        }
      }
      if (n == 1){ n = 2;}
      stddev = (stddev/div)*n/(n-1);
      htot1->SetBinContent(nrebins+1, sum/div);

      Double_t err = sqrt(abs(sum)/div);
      htot1->SetBinError(nrebins+1,sqrt(stddev + err*err));
      if(sum!=0 && htot1->GetBinCenter(nrebins+1)>threshold/3.8){
        _h1tot[0].push_back(sum/div);
        _h1tot[1].push_back(htot1->GetBinError(nrebins+1));
        _hpe[0].push_back(htot1->GetBinCenter(nrebins+1));
        _hpe[1].push_back(htot1->GetBinWidth(nrebins+1)/2.);
      }

      nrebins++;

      sum = 0;
      div = 0;
      stddev = 0;
      n = 0;
      ndiv = 0;
    }


  }
}

void configHisto(TH1D *h, Color_t color){
  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.7);
  h->SetLineWidth(2);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->GetYaxis()->SetTitle("TOT (ns)");
  h->GetYaxis()->SetRangeUser(10,800);
  h->GetXaxis()->SetTitle("Photo-electrons");
  h->GetXaxis()->SetRangeUser(0.5*threshold/3.4,3000);
}

void configHisto2(TH2D *h){
  h->GetXaxis()->SetTitle("TOT (ns)");
  h->GetXaxis()->SetRangeUser(10,800);
  h->GetYaxis()->SetTitle("Photo-electrons");
  h->GetYaxis()->SetRangeUser(0.5*threshold/4.,3000);
  h->SetStats(kFALSE);
}


void configGraphs(TGraph *g, Color_t c){
  
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.7);
  g->SetMarkerColor(c);
  g->SetLineColor(c);
  g->SetLineWidth(2);
  g->GetXaxis()->SetTitle("TOT (ns)");
  g->GetXaxis()->SetRangeUser(10,800);
  g->GetYaxis()->SetTitle("Photo-electrons");
  g->GetYaxis()->SetRangeUser(0.5*threshold/4.,3000);
}
void reconstruct(TF1 *fcorr){
  Int_t n = gtot_corr.size();
  gpe_corr.resize(n);
  gpe_corr_r.resize(n);
  er_gpe_corr.resize(n);
  for (Int_t i = 0; i < n; i++) {
    gpe_corr[i] = fcorr->Eval(r3->Gaus(gtot_corr[i],1));
    er_gpe_corr[i] = sqrt(gpe_corr[i]);

    hdev_sat->Fill((gpe_sat[i] - gpe[i])/gpe[i]);
    gpe_corr_r[i] = (gpe_corr[i] - gpe[i])/gpe[i];
    hdev->Fill((gpe_corr[i] - gpe[i])/gpe[i]);
    // htot_sat->Fill(gpe_corr[i], gtot_sat[i]);
  }
  TCanvas *cc = new TCanvas("cc","cc",1920,0,1920,1080);
  cc->cd();
  TPad *p1 = new TPad("","", 0, 0, 1, 1);
  TPad *p2 = new TPad("","", 0, 0, 1, 1);
  Double_t ysplit = 0.3;
  p1->SetBottomMargin(ysplit);
  p2->SetTopMargin(1-ysplit);
  p1->SetFillStyle(0);
  p2->SetFillStyle(0);


  cc->cd(); p2->Draw(); p2->cd();
  TGraph *gcr = new TGraph(n, &gpe[0], &gpe_corr_r[0]);
  TGraph *gsr = new TGraph(n, &gpe[0], &gpe_sat_r[0]);
  configGraphs(gcr, kBlack);
  configGraphs(gsr, kRed);
  gcr->SetMarkerColorAlpha(kBlack,0.01);
  gcr->GetXaxis()->SetTitle("Photo-electrons (origial)");
  gcr->GetXaxis()->SetRangeUser(500,2400);
  gcr->GetYaxis()->SetTitle("Ratio");
  gcr->GetYaxis()->SetRangeUser(-1,1);
  gcr->GetYaxis()->SetNdivisions(205);

  gcr->Draw("AP");
  gsr->Draw("SAME P");

  p2->SetGrid(1,1);
  p1->cd();


  cc->cd(); p1->Draw(); p1->cd();
  TGraphErrors *gbad = new TGraphErrors(n, &gpe[0], &gpe_sat[0], &er_gpe[0], &er_gpe_sat[0]);//, &er_gpe[0], &er_gpe_sat[0]);
  TGraphErrors *gcorr = new TGraphErrors(n, &gpe[0], &gpe_corr[0], &er_gpe[0], &er_gpe_corr[0]);//, &er_gpe[0], &er_gpe_corr[0]);




  
  configGraphs(gbad,kRed);
  configGraphs(gcorr,kBlack);

  gcorr->SetMarkerColorAlpha(kBlack,0.01);
  gbad->GetXaxis()->SetTitle("Photo-electrons (origial)");
  gbad->GetXaxis()->SetRangeUser(500,2400);
  gbad->GetXaxis()->SetTitleSize(0);
  gbad->GetXaxis()->SetLabelSize(0);
  gbad->GetYaxis()->SetTitle("Photo-electrons");
  gbad->GetYaxis()->SetRangeUser(500,3000);
  gbad->SetTitle("Saturated events");
  gcorr->SetTitle("TOT Corrected events");

  gbad->Draw("AP X");
  gcorr->Draw("SAME P X");
  TF1 *fx = new TF1("fx","[0]*x",0,2400);
  fx->FixParameter(0,1);
  fx->SetTitle("y = x");
  fx->SetRange(0,3000);
  gcorr->Fit("fx","R");
  TLegend *lg = new TLegend(0.1,0.7,0.5,0.9);
  lg->AddEntry(gbad);
  lg->AddEntry(gcorr);
  lg->Draw();



  TCanvas *chdev = new TCanvas("chdev","chdev",1920,0,1920,1080);
  hdev_sat->SetLineColor(kRed);
  hdev_sat->SetLineWidth(2);
  hdev->SetLineWidth(2);
  hdev_sat->GetYaxis()->SetTitle("# of events");
  hdev_sat->GetXaxis()->SetTitle("Relative deviation");
  hdev_sat->SetTitle("#frac{Saturated - original}{original}");
  hdev->SetTitle("#frac{Corrected - original}{original}");
  hdev_sat->Draw("");
  hdev->Draw("SAMES");

  chdev->BuildLegend();
  // cc->Print("graphs/TOT_graphs.root");
  chdev->Print(Form("graphs/TOT_relative_amp_th%d.root", (int)threshold));



}
void tot_test_amp(){
  gStyle->SetOptTitle(0);
  if(threshold == 200) minFit = 200;
  if(threshold == 300) minFit = 120;
  if(threshold == 400) minFit = 100;
  if(threshold == 600) minFit = 80;
  if(threshold == 800) minFit = 80;
  if(threshold == 1000) minFit = 100;
  if(mychannel == "Ch3"){minInt = 10280; maxInt = 10700;}
  for(Int_t i = 0; i<n; i++){
    files[i] = files[i]+conc;
    z[i] = new ANALYZER(Form("z%.2f",volts[i]));
    z[i]->setAnalyzer(files[i].c_str());
    z[i]->setChannel(mychannel.c_str());
    Int_t kch = z[i]->kch;
    saturation_level = saturations[z[i]->getIdx()];
    sphe = sphes[z[i]->getIdx()];

    Int_t nentries = (just_a_test) ? 100 : z[i]->nentries;
    for(Int_t j = 0; j < nentries; j++){
      z[i]->getWaveform(j,kch);
      z[i]->applyDenoise(filter);
      if(z[i]->ch[kch].selection!=0) continue;
      // z[i]->applyDenoise(16);
      Double_t tot = z[i]->getTOT(kch, minInt, 15000, threshold, 12);
      z[i]->integrate(kch, minInt-200, minInt);
      if (z[i]->temp_max > 20 || tot <= 12){
        continue;
      }
      z[i]->integrate(kch,minInt,maxInt);

      Double_t pe = z[i]->temp_max/sphe;
      htot->Fill(pe, tot);
      hQ->Fill(tot, pe);
      if (z[i]->temp_max >= saturation_level){
        gpe.push_back(pe);
        er_gpe.push_back(sqrt(pe));
        gtot.push_back(tot);
        er_gtot.push_back(8);

        z[i]->getWaveform(j,kch);
        for (Int_t k = 0; k < memorydepth; k++) {
          z[i]->ch[kch].wvf[k] = (z[i]->ch[kch].wvf[k] >= saturation_level) ? r3->Gaus(saturation_level,0.5*sqrt(saturation_level)) : z[i]->ch[kch].wvf[k];
        }
        z[i]->applyDenoise(filter);
        tot = z[i]->getTOT(kch, minInt, 15000, threshold, 12);
        z[i]->integrate(kch,minInt,maxInt);
        Double_t prev_pe = pe;
        pe = z[i]->temp_max/sphe;
        gpe_sat.push_back(pe);
        gpe_sat_r.push_back((pe - prev_pe)/prev_pe);
        er_gpe_sat.push_back(sqrt(pe));
        // gtot_sat.push_back(tot);
        er_gtot_sat.push_back(8);
        gtot_corr.push_back(tot);
        er_gtot_corr.push_back(8);
      }
      else{
        htot_before_sat->Fill(pe, tot);
        hQ_before_sat->Fill(tot, pe);
      }
      htot_sat->Fill(pe, tot);
      hQ_sat->Fill(tot, pe);
      hshow_amp_vs_pe->Fill(pe,z[i]->temp_max);
    }

  }

  setup_histograms(htot, htot1, h1tot, hpe);
  setup_histograms(htot_sat,htot1_sat, h1tot_sat, hpe_sat);
  setup_histograms(htot_before_sat,htot1_before_sat, h1tot_before_sat, hpe_before_sat);


  TGraphErrors *gfull = new TGraphErrors(hpe[0].size(), &h1tot[0][0], &hpe[0][0], &h1tot[1][0], &hpe[1][0]);
  TGraphErrors *g_sat = new TGraphErrors(hpe_sat[0].size(), &h1tot_sat[0][0], &hpe_sat[0][0], &h1tot_sat[1][0], &hpe_sat[1][0]);
  TGraphErrors *g_before_sat = new TGraphErrors(hpe_before_sat[0].size(), &h1tot_before_sat[0][0], &hpe_before_sat[0][0], &h1tot_before_sat[1][0], &hpe_before_sat[1][0]);
  configGraphs(gfull, kBlack);
  configGraphs(g_sat, kRed);
  configGraphs(g_before_sat, kBlack);

  TF1 *f = new TF1("f","expo(0)",0,2500);
  TF1 *fexpo0 = new TF1("fexpo0","expo(0)",0,2500);
  // Double_t param[4] = {3.3, 0.007, 3.4, -0.004};
  Double_t param[4] = {2.8, 0.007, 3.4, -0.006};


  f->SetParameters(param[0],param[1],param[2],param[3]);
  fexpo0->SetParameters(param[0],param[1]);

  TCanvas *c1 = new TCanvas("c1","c1",1920,0,1920,1080);
  configHisto2(hQ);
  hQ->Draw("colz");
  gfull->Draw("SAME ZP");
  gfull->Fit("fexpo0","W0", "", minFit,3000);
  gfull->Fit("fexpo0","0", "", minFit,3000);
  gfull->Fit("fexpo0","0", "", minFit,3000);
  f->SetParameters(fexpo0->GetParameter(0), fexpo0->GetParameter(1));
  // gfull->Fit("f","W");
  gfull->Fit("f","");
  TF1 *fperf = (TF1*)f->Clone("fperf");
  for(Int_t l = 0; l < 4; l++){
    cout << fperf->GetParameter(l) << " ";
  }
  cout << "\n";
  fperf->SetLineColor(kGreen);
  fperf->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2","c2",1920,0,1920,1080);
  configHisto2(hQ_sat);
  hQ_sat->Draw("colz");
  g_sat->Draw("SAME ZP");

  TCanvas *c3 = new TCanvas("c3","c3",1920,0,1920,1080);
  hshow_amp_vs_pe->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4",1920,0,1920,1080);
  configHisto2(hQ_before_sat);
  hQ_before_sat->Draw("colz");
  g_before_sat->Draw("SAME ZP");
  fexpo0->SetParameters(param[0],param[1]);
  g_before_sat->Fit("fexpo0","W", "", minFit, 3000);
  g_before_sat->Fit("fexpo0","", "", minFit, 3000);
  f->SetParameters(param[0],param[1],param[2],param[3]);
  f->SetParameters(fexpo0->GetParameter(0), fexpo0->GetParameter(1));
  g_before_sat->Fit("f","W", "", minFit, 3000);
  g_before_sat->Fit("f", "", "", minFit, 3000);
  g_before_sat->Fit("f", "", "", minFit, 3000);
  f->Draw("SAME");
  fperf->Draw("SAME");
  TLegend *l4 = new TLegend(0.1,0.7,0.5,0.9);
  l4->AddEntry(f,"Fit: exp([0]+[1]*x) + exp([2]+[3]*x)", "l");
  l4->AddEntry(fperf,"Fit without saturation", "l");
  l4->Draw();


  // c1->Print("graphs/TOT_full.root");
  // c2->Print("graphs/TOT_sat.root");
  // c4->Print("graphs/TOT_before_sat.root");
  reconstruct(f);


}
