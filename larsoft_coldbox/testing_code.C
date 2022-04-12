// analysistree->cd();anatree->GetMaximum("no_hits_stored");anatree->GetMaximum("ntracks_pandoraTrack")
#define maxHits 1617
#define maxTracks 14
#include "class/analyzer.C"
void testing_code(){
  
  DATA data("ana_hist_429_1.root");
  DATA data2("ana_hist_429_2.root");
  DATA data3("ana_hist_429_3.root");
  // DATA data4("ana_hist_455_1.root");
  // data.read_tree(); // dont need this anymore


  data.selection_phi();
  data.FillHistoWithData(data);
  data.FillHistoWithData(data2);
  data.FillHistoWithData(data3);
  // data.FillHistoWithData(data4);
  data.plot_tracks();
  
  data.test_some_events();

  TF1 *fmu = new TF1("fmu","[0]*pow(cos(x),[1])*sin(x)+[2]",90,180);
  fmu->SetParameters(1,2,0.1);
  fmu->SetParName(0,"A");
  fmu->SetParName(1,"n");
  fmu->SetParName(2,"B");

  // for(Int_t i = 1; i<=110; i++)
  //   {
  //     data.test_some_events(i,2,0);
  //   }
  



}
