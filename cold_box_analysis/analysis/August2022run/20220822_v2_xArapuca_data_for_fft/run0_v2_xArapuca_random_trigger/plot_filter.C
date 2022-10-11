#define memorydepth 125000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void plot_filter(Int_t myevent = 2){
  TFile *f = new TFile("analyzed.root","READ");
  TTree *t1 = (TTree*)f->Get("t1");

  vector<Int_t> channels = {1};
  Int_t n = channels.size();
  ADC_DATA ch;
  TBranch *bch = t1->GetBranch("Ch1");
  bch->SetAddress(&ch);

  TCanvas* c = new TCanvas();

  Double_t dtime = 4;
  Double_t freq = 250;
  TH1D *hsignal = new TH1D("hsignal", "hsignal", memorydepth, 0, dtime * memorydepth);

  WIENER wo("wo",dtime,freq,1e-9,1e6,memorydepth);
  WIENER ws("ws",dtime,freq,1e-9,1e6,memorydepth);


  bch->GetEvent(myevent);
  for (Int_t j = 0; j < memorydepth; j++) {
    hsignal->SetBinContent(j+1,ch.wvf[j]);
  }
  ws.fft(hsignal);
  wo.fft(hsignal);
  ws.apply_filter(20);
  ws.backfft(ws);

  hsignal->SetLineColor(kRed);
  hsignal->Draw();
  ws.hwvf->Draw("SAME");


  TCanvas *c2 = new TCanvas();
  wo.hfft->SetLineColor(kRed);
  wo.hfft->Draw();
  ws.hfft->Draw("SAME");
}
