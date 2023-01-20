#define memorydepth 250000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

vector<string> devices = {"0","0","miniArapuca 37V (A1ch1)","miniArapuca 47V (Argon4)", "xArapuca v4 ch1 (DCemArgon4)", "xArapuca v4 ch2 (DCemArgon4)", "xArapuca v5 ch1 (DCemSimp3)", "xArapuca v5 ch2 (DCemSimp3)"};
vector<Double_t> saturations = {0, 0, 13000, 14000, 12000, 12000, 12000, 100};

ANALYZER zall("zall");
ANALYZER znov5("znov5");
ANALYZER znopof("znopof");

void getffts(ANALYZER &z, THStack *hs){
  for(Int_t k = 0; k < z.nchannels; k++){
  // for(Int_t k = 0; k < 1; k++){
    cout << "Channel " << k << endl;
    z.kch = k;
    if(z.schannel[k] == "Ch7") continue;
    z.getSelection(""); // need this to reset selection from previous channel
    z.selectByAmplitude(16,0,z.dtime*memorydepth, saturations[z.channels[z.kch]]);
    cout << "done.." << endl;
    z.averageFFT(0, "use_mine", true, 16);
    z.hfft[k]->SetTitle(devices[z.channels[k]].c_str());
    hs->Add(z.hfft[k]);
  }
}
void ffts(){
  zall.setAnalyzer("run19_all_devices_random_trigger/analyzed.root");
  znov5.setAnalyzer("run20_v5_turned_off_random_trigger/analyzed.root");
  znopof.setAnalyzer("run21_v5_and_v4_turned_off_random_trigger/analyzed.root");

  THStack *hs = new THStack("hs","hs");
  THStack *hs2 = new THStack("hs2","hs2");
  THStack *hs3 = new THStack("hs3","hs3");

  getffts(zall,hs);
  getffts(znov5,hs2);
  getffts(znopof,hs3);

  // TCanvas *c1 = new TCanvas("c1","c1");
  // hs->Draw("nostack PLC");

  // TCanvas *c2 = new TCanvas("c2","c2");
  // hs2->Draw("nostack PLC");

  // TCanvas *c3 = new TCanvas("c3","c3");
  // hs3->Draw("nostack PLC");

  TCanvas *c4 = new TCanvas("c4", "c4");
  Int_t a = 0;
  zall.hfft[a]->Draw("");
  zall.hfft[a]->SetLineWidth(2);
  zall.hfft[a]->SetLineColor(kBlack);
  znov5.hfft[a]->Draw("SAMES");
  znov5.hfft[a]->SetLineWidth(2);
  znov5.hfft[a]->SetLineColor(kBlue);
  znopof.hfft[a]->Draw("SAMES");
  znopof.hfft[a]->SetLineWidth(2);
  znopof.hfft[a]->SetLineColor(kRed);

  
}
