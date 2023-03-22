#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void time_res(){
  ANALYZER z("z");
  z.setAnalyzer();
  TH1D *hres = new TH1D("hres","hres",500,-100,100);
  vector<Double_t> vmax3;
  vector<Double_t> vmax4;
  for(Int_t i = 0; i < z.nentries; i++){
    z.kch = 3;
    z.getWaveform(i);
    z.applyDenoise(16);
    Double_t max3 = z.getMaximum(10250, 10450);
    if(max3 < 20) continue;
    z.triggerTime({0,0}, {10250, 10450},10);
    Double_t pos3 = z.temp_pos;
    z.kch = 4;
    z.getWaveform(i);
    z.applyDenoise(16);
    // if(z.ch[z.kch].selection!=0) continue;
    Double_t max4 = z.getMaximum(10250, 10450);
    if(max4 >= 13000) continue;
    if(max4 < 100) continue;
    z.triggerTime({0,0}, {10250, 10450},240);
    Double_t pos4 = z.temp_pos;
    Double_t difftime = pos3 - pos4;
    if(pos4 == 0 || pos3 == 0) continue;
    hres->Fill(difftime);

    vmax3.push_back(max3);
    vmax4.push_back(max4);
    if(difftime < -20){
      cout << i << " " <<  difftime << " " << pos3 << " " << pos4 << endl;
    }
  }

  hres->Draw();
  hres->GetXaxis()->SetTitle("Time difference (ns)");
  hres->GetYaxis()->SetTitle("# of events");

  TCanvas *c2 = new TCanvas("c2", "c2",1920,0,1920,1080);
  TGraph *g = new TGraph(vmax3.size(),&vmax3[0],&vmax4[0]);
  g->Draw("AP*");

}
