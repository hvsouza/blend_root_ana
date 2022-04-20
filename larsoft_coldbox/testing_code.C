// analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
#define maxHits 1617
#define maxTracks 14
#include "class/analyzer.C"
void testing_code(){
  
  DATA data("ana_hist_429_1.root");
  DATA data2("ana_hist_429_2.root");
  DATA data3("ana_hist_429_3.root");
  DATA data4("ana_hist_429_4.root");
  DATA data5("ana_hist_429_5.root");
  DATA data6("ana_hist_429_6.root");
  // DATA dataextra("ana_hist_455_1.root");
  // data.read_tree(); // dont need this anymore


  // data.selection_phi();
  data.FillHistoWithData(data);
  data.FillHistoWithData(data2);
  data.FillHistoWithData(data3);
  data.FillHistoWithData(data4);
  data.FillHistoWithData(data5);
  data.FillHistoWithData(data6);
  // data.FillHistoWithData(dataextra);
  data.plot_tracks();
  
  data.test_some_events();

  TF1 *fmu = new TF1("fmu","([0]*pow(cos(x*TMath::Pi()/180),[1])*sin(x*TMath::Pi()/180))*cos(x*TMath::Pi()/180)",0,90);
  fmu->SetParameters(1,2// ,0.1,0.1
                     );
  fmu->SetParName(0,"A");
  fmu->SetParName(1,"n");
  data.htheta->Fit("fmu","R");
  TRandom3 *rd = new TRandom3(93);
  TH1D *hres = new TH1D("hres","hres",data.nbins_htheta,0,90);
  // hres->FillRandom("fmu",data.htheta->GetEntries()*5,rd);
  hres->FillRandom("fmu",1000000,rd);
  hres->Scale(1./hres->Integral("width"));
  hres->SetLineWidth(2);
  hres->SetLineColor(kRed);
  // hres->Draw("SAME");
  fmu->Draw("SAME");


  // for(Int_t i = 1; i<=110; i++)
  //   {
  //     data.test_some_events(i,2,0);
  //   }
  



}
