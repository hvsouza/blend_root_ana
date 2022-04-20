// analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
#define maxHits 11125
#define maxTracks 62
#include "/home/henrique/Dropbox/APC_Paris/Root/larsoft_coldbox/class/analyzer.C"
void testing_code(){
  
  DATA data("429/429_cb_reco.root");
  // data.maxPlotEvents = 500;
  // data.selection_phi();
  data.FillHistoWithData(data);
  data.plot_tracks();
  
  data.test_some_events();

  TF1 *fmu = new TF1("fmu","([0]*pow(cos(x*TMath::Pi()/180),[1])*sin(x*TMath::Pi()/180))*cos(x*TMath::Pi()/180)",0,90);
  fmu->SetParameters(1,2// ,0.1,0.1
                     );
  fmu->SetParName(0,"A");
  fmu->SetParName(1,"n");
  data.htheta->Fit("fmu","RI");
  // TRandom3 *rd = new TRandom3(93);
  // TH1D *hres = new TH1D("hres","hres",data.nbins_htheta,0,90);
  // // hres->FillRandom("fmu",data.htheta->GetEntries()*5,rd);
  // hres->FillRandom("fmu",1000000,rd);
  // hres->Scale(1./hres->Integral("width"));
  // hres->SetLineWidth(2);
  // hres->SetLineColor(kRed);
  // hres->Draw("SAME");
  fmu->Draw("SAME");


  // for(Int_t i = 1; i<=110; i++)
  //   {
  //     data.test_some_events(i,2,0);
  //   }
  



}
