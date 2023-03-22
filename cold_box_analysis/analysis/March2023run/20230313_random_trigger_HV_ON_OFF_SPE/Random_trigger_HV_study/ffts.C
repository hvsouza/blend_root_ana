#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

vector<Int_t> channels = {0,3,4,5,6};
vector<string> devices = {"miniArapuca 37V (A1ch1)", "xArapuca v4 ch2 (DCem1.2)", "xArapuca v4 ch1 (DCem1.2)", "xArapuca v5 ch2 (DCemArgon4)", "xArapuca v5 ch1 (DCemArgon4)"};
vector<Double_t> saturations = {12000, 12000, 12000, 12000, 12000, 12000};

vector<string> description = {"HV ON, CRP bias on", "HV ON, CRP bias off", "HV OFF", "HV OFF (no Argon2x2)", "HV OFF - HD electronic off", "HV ON again, CRP bias on", "HV ON again, power cycle", "HV ON, v4 working"};
vector<Color_t> mycolor = {kRed, kViolet, kBlack, kBlue, kAzure+1, kGreen, kGreen+3, kYellow+2};

vector<string> files = {
                        "./HV_ON_start/run1_argon2x2_and_v4/analyzed.root",
                        "./HV_ON_start/run2_argon2x2_and_v4_and_v5/analyzed.root"//,
                        "./HV_ON_start/run3_argon2x2_and_v4_and_v5_CRP_bias_off/analyzed.root",
                        "./HV_OFF/run4_argon2x2_and_v4_and_v5_CRP_bias_off_HV_OFF/analyzed.root",
                        "./HV_OFF/run5_v4_and_v5_CRP_bias_off_HV_OFF/analyzed.root",
                        "./HV_OFF/run6_argon2x2_v4_and_v5_CRP_bias_off_and_electronics_off_HV_OFF/analyzed.root",
                        "./HV_ON_after/run7_argon2x2_v4_and_v5_after_HV_ON/analyzed.root",
                        "./HV_ON_after/run8_argon2x2_v4_and_v5_power_cycle_HD_system_on_HV_ON/analyzed.root",
                        "./HV_ON_after/run9_argon2x2_v4_and_v5_power_increased_HD_system_on_HV_ON/analyzed.root"
                        };

vector<ANALYZER *> z;
Int_t nfiles;

void getffts(ANALYZER &z){
  Double_t factor = pow(2,14);
  Int_t maxev = 0;
  for(Int_t k = 0; k < (int)channels.size(); k++){
    z.hfft[k] = nullptr;
  // for(Int_t k = 0; k < 1; k++){
    if (k > 2){ continue; }
    if(!z.setChannel(Form("Ch%d.", channels[k]))) continue;

    cout << z.kch << ": " << z.schannel[k] << endl;
    if(z.schannel[k] == "Ch7.") continue;
    z.getSelection(""); // need this to reset selection from previous channel
    // z.selectByAmplitude(16,0,z.dtime*z.n_points, saturations[k]);
    cout << "done selection: " << z.lev->GetN() << endl;
    z.averageFFT(maxev, "", true, 0, 16);
    // z.h->Rebin(4);
    // z.hfft[k]->SetTitle(devices[k].c_str());
  }
}
void plot_some_graph(TCanvas *c, Int_t a = 0){

  bool draw = false;
  for(Int_t i = 0; i < nfiles; i++){
    if(!z[i]->hfft[a]) continue;
    if(!draw){
      draw = true;
      z[i]->hfft[a]->Draw("");
    }
    else z[i]->hfft[a]->Draw("SAME");
    z[i]->hfft[a]->SetLineWidth(3);
    z[i]->hfft[a]->SetLineColor(mycolor[i]);
    z[i]->hfft[a]->SetTitle(Form("%s - %s",devices[a].c_str(), description[i].c_str()));
  }
  c->BuildLegend();
}
void ffts(){
  nfiles = files.size();
  z.resize(nfiles);
  for (int i = 0; i < nfiles; i++) {
    z[i] = new ANALYZER(Form("z%d",i));
    z[i]->setAnalyzer(files[i]);
    cout << files[i] << endl;
    getffts(*z[i]);
    cout << "\n";
  }

  TCanvas *c4 = new TCanvas("c4", "c4");
  plot_some_graph(c4,0);

  TCanvas *c5 = new TCanvas("c5", "c5");
  plot_some_graph(c5,1);

  TCanvas *c6 = new TCanvas("c6", "c6");
  plot_some_graph(c6,2);

  // TCanvas *c7 = new TCanvas("c7", "c7");
  // plot_some_graph(c7,3);

  // TCanvas *c8 = new TCanvas("c8", "c8");
  // plot_some_graph(c8,4);
}
