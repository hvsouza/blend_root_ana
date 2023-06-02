#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

WIENER *wspe = nullptr;


void getSPEFFT(){

  cout << "Computing SPE..." << endl;

  ANALYZER zpe("zpe");
  zpe.setAnalyzer("./sphe_waveforms_Ch6.root");
  // TGraph *gsphe = (TGraph*)zpe.f->Get("mean");
  // zpe.getWaveFromGraph(gsphe);
  //
  vector<Double_t> avgspe(zpe.n_points);
  Int_t navg = 0;
  zpe.getSelection("Ch6.selection==1");
  cout << "Got selected by fit" << endl;
  zpe.selectByAmplitude(16, 0, 5200, 25, "higher"); // remove pretrigger with pulse
  zpe.selectByAmplitude(16, 5450, 10000, 25, "higher"); // remove pretrigger with pulse
  cout << "Applied cut for amplitude" << endl;
  zpe.selectByAmplitude(16, 5500, 10000, -25, "lower"); // remove some undershoot
  cout << "Applied cut for undershoot" << endl;

  for (int i = 0; i < zpe.lev->GetN(); i++) {
    Int_t iev = zpe.lev->GetEntry(i);
    zpe.getWaveform(iev);
    for(Int_t j = 0; j < zpe.n_points; j++){
      avgspe[j] += zpe.ch[0]->wvf[j];
    }
    navg += 1;

  }
  cout << "Valid SPEs = " << navg << endl;
  for(Int_t i = 0; i < zpe.n_points; i++){
    wspe->hwvf->SetBinContent(i+1, avgspe[i]/navg);
  }
  wspe->fft(wspe->hwvf);
  zpe.persistence_plot(500,-100,500,16,"use_mine");
  wspe->hwvf->Draw("SAME");

}
void deconvolve(){


  Double_t spe_peak = 29.6;
  vector<Double_t> range_av = {200,300};


  ANALYZER z("z");
  z.setAnalyzer("./analyzed.root");
  z.setChannel("Ch6.");

  Int_t n_points = z.n_points;
  Int_t kch = z.kch;
  Double_t *wvf = nullptr;

  vector<Double_t> avg_wvf(n_points, 0);

  wspe = new WIENER("wspe", 4, 250, 1e-9, 1e6, n_points);
  getSPEFFT();
  Int_t nav = 0;
  for(Int_t i = 0; i < z.nentries; i++){
  // for(Int_t i = 0; i < 100; i++){
    if(i%200==0) cout << "computing event " << i << " of " << z.nentries << "\r" << flush;
    z.getWaveform(i);
    wvf = z.ch[kch]->wvf;
    if(z.ch[kch]->selection == 1) continue;
    Double_t max = z.getMaximum(2600, 3000);
    if(max/spe_peak < range_av[0] || max/spe_peak > range_av[1]) continue;
    z.getFFT();
    z.drawGraph("");
    z.w->deconvolve(*z.w, *wspe, 5, "gaus");
    // z.w->hwvf->SetLineColor(kRed);
    // z.w->hwvf->Draw("SAME");
    // break;
    for(Int_t j = 0; j < n_points; j++){
      // avg_wvf[j] +=  wvf[j];
      avg_wvf[j] +=  z.w->hwvf->GetBinContent(j+1);
    }
    nav++;
  }
  Double_t *time = new Double_t[n_points];
  for(Int_t j = 0; j < n_points; j++){
    avg_wvf[j] = avg_wvf[j]/nav;
    time[j] = j*z.dtime;
  }
  cout << "\n";
  cout << "Found " << nav << " waveforms" << endl;
  TGraph *gwvf = new TGraph(n_points, time, &avg_wvf[0]);
  gwvf->Draw("AL");
  gwvf->GetXaxis()->SetTitle("Time (ns)");
  gwvf->GetYaxis()->SetTitle("Amplitude (A.U.)");


}
