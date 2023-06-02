#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

WIENER *wspe = nullptr;
WIENER *wnoise = nullptr;
WIENER *wexpect = nullptr;
WIENER *wg = nullptr;
WIENER *wsignal = nullptr;

void getExpected(Double_t amp){

  TF1 *fexp = new TF1("fexp", "[0]*exp(-x/[1]) + [2]*exp(-x/[3])", 0, 10000);
  //assuming that at t = 0, the amplitude is only due to Af:
  // Af = amp;
  // Af*tau_f = amp*7;
  // As*tau_s = 3 * amp*7 (75% slow)
  // As = 3*amp*7 / tau_s
  fexp->SetParameters(amp, 7, 3*amp*7/1600, 1600);

  for(Int_t i = 0; i < wexpect->npts; i++){
    wexpect->hwvf->SetBinContent(i+1, fexp->Eval(i*4));
  }

  wexpect->fft(wexpect->hwvf);

  // wexpect->hfft->Draw();

}

void getNoise(){

  ANALYZER zn("zn");
  zn.setAnalyzer("./sphe_waveforms_Ch6.root");
  vector<Double_t> avg_noise(zn.n_points,0);
  int aux = 0;
  int navg = 0;
  int nnavg = 0;
  Double_t *dwvf = new Double_t[zn.n_points];
  WIENER *wtemp = new WIENER("wtemp", 4, 250, 1e-9, 1e6, zn.n_points);
  Int_t groupn = 200;
  for(Int_t i = 0; i < zn.nentries; i++){
    zn.getWaveform(i);
    zn.applyDenoise(16, zn.ch[0]->wvf, dwvf);
    bool good_noise = true;;
    for (int j = 0; j < zn.n_points; j++) {
      if((dwvf[j]) >= 18 || dwvf[j] < -35){
        good_noise = false;
        break;
      }
    }
    if(!good_noise) continue;
    zn.lev->Enter(i);
    for (int j = 0; j < zn.n_points; j++) {
      avg_noise[j] += zn.ch[0]->wvf[j];
    }
    navg++;
    if(navg>groupn || i == zn.nentries - 1){
      nnavg++;
      for (int j = 0; j < zn.n_points; j++) {
        wtemp->hwvf->SetBinContent(j+1, avg_noise[j]/navg);
      }
      wtemp->fft(wtemp->hwvf);
      avg_noise.clear();
      navg = 0;
      wnoise->add_fft(*wtemp);
    }
  }
  wnoise->scale(1./nnavg);
  wnoise->recompute_hist();
  // wnoise->hPSD->Draw("hist");
  // zn.persistence_plot(500,-100,500,16,"use_mine");

  cout << "Noise waveforms found = " << nnavg*groupn << endl;
}
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
void deconvolve_wiener(){

  THStack *hs = new THStack("hs","hs");
  TLegend *lg = new TLegend(0.8,0.8,1,1);

  Double_t spe_peak = 29.6;
  vector<Double_t> range_av = {200,300};


  ANALYZER z("z");
  z.setAnalyzer("./analyzed.root");
  z.setChannel("Ch6.");

  Int_t n_points = z.n_points;
  Int_t kch = z.kch;
  Double_t *wvf = z.ch[kch]->wvf;

  vector<Double_t> avg_wvf(n_points, 0);

  wspe = new WIENER("wspe", 4, 250, 1e-9, 1e6, n_points);
  getSPEFFT();
  hs->Add(wspe->hPSD,"hist");
  lg->AddEntry(wspe->hPSD,"SPE");

  wnoise = new WIENER("wnoise", 4, 250, 1e-9, 1e6, n_points);
  getNoise();
  hs->Add(wnoise->hPSD,"hist");
  lg->AddEntry(wnoise->hPSD,"Noise");

  double amp = (range_av[1] + range_av[0])/2.;
  wexpect = new WIENER("wexpect", 4, 250, 1e-9, 1e6, n_points);
  getExpected(amp);
  hs->Add(wexpect->hPSD,"hist");
  lg->AddEntry(wexpect->hPSD, "LAr response");

  wg = new WIENER("wg", 4, 250, 1e-9, 1e6, n_points);
  wg->wienerGenFilter(*wspe, *wnoise, *wexpect, 2);
  hs->Add(wg->hPSD,"hist");
  lg->AddEntry(wg->hPSD,"Wiener filter");
  // wg->hfft->Draw();

  wsignal = new WIENER("wsignal", 4, 250, 1e-9, 1e6, n_points);

  Int_t nav = 0;
  for(Int_t i = 0; i < z.nentries; i++){
  // for(Int_t i = 0; i < 100; i++){
    if(i%200==0) cout << "computing event " << i << " of " << z.nentries << "\r" << flush;
    z.getWaveform(i);
    if(z.ch[kch]->selection == 1) continue;
    Double_t max = z.getMaximum(2600, 3500);
    if(max/spe_peak < range_av[0] || max/spe_peak > range_av[1]) continue;
    for(Int_t j = 0; j < n_points; j++){
      avg_wvf[j] +=  z.ch[kch]->wvf[j];
    }
    nav++;
  }
  Double_t *time = new Double_t[n_points];
  for(Int_t j = 0; j < n_points; j++){
    avg_wvf[j] = avg_wvf[j]/nav;
    wsignal->hwvf->SetBinContent(j+1, avg_wvf[j]);
    time[j] = j*z.dtime;
  }
  cout << "\n";
  cout << "Found " << nav << " waveforms" << endl;


  // TGraph *gwvf = new TGraph(n_points, time, &avg_wvf[0]);
  // gwvf->Draw("AL");
  // gwvf->GetXaxis()->SetTitle("Time (ns)");
  // gwvf->GetYaxis()->SetTitle("Amplitude (A.U.)");


  TCanvas *c3 = new TCanvas();

  wsignal->fft(wsignal->hwvf);
  TH1D *hffts = (TH1D*)wsignal->hPSD->Clone("hffts");
  hs->Add(hffts, "hist");
  lg->AddEntry(hffts, "Average signal");
  wsignal->fft(wsignal->hwvf);
  wsignal->applyWienerFilter(*wg);
  wsignal->hwvf->Draw();
  hs->Add(wsignal->hPSD, "hist");
  lg->AddEntry(wsignal->hPSD, "Deconvoluted");

  TCanvas *c2 = new TCanvas();
  hs->Draw("nostack plc");
  hs->GetXaxis()->SetTitle("Frequency (MHz)");
  hs->GetYaxis()->SetTitle("PSD (ADCs^{2} MHz^{-1})");
  lg->Draw();

}
