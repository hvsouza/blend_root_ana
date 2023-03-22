#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void write_timestemp(){
  ANALYZER z("z");
  z.setAnalyzer();

  ofstream fout;
  fout.open("PDS_TriggerTimes.txt",ios::out);

  fout << "# event,\ttime (s)\n";
  // for(Int_t i = 0; i < 265; i++){
  for(Int_t i = 0; i < z.nentries; i++){
    z.getWaveform(i,0);
    fout << (int)z.ch[0]->event << ",\t" << z.ch[0]->time << "\n";
  }
  fout.close();
}
