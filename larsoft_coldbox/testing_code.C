#include "class/analyzer.C"
void testing_code(){
  DATA data;

  data.read_tree();
  data.plot_tracks();
  
  data.test_some_events();
  
  // for(Int_t i = 1; i<=110; i++)
  //   {
  //     data.test_some_events(i,2,0);
  //   }
  



}
