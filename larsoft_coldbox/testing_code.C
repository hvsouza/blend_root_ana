#include "class/analyzer.C"
void testing_code(){
  DATA data("ana_hist.root");

  // data.read_tree(); // dont need this anymore
  data.selection_phi();
  data.plot_tracks();
  
  data.test_some_events();
  
  // for(Int_t i = 1; i<=110; i++)
  //   {
  //     data.test_some_events(i,2,0);
  //   }
  



}
