// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer.C
// Created: 2022
// ________________________________________ //
//
#include "MYCODES.h"



class ANALYZER{
  public:

    TFile *f = nullptr;
    TTree *t1 = nullptr;
    Int_t nentries = 0;
    vector<TBranch*> b;
    vector<ADC_DATA> ch;
    Int_t nchannels = 0;
    vector<Int_t> channels = {1,2};
    vector<string> schannel;
    Double_t dtime = 4;
    string myname;
    string filename = "";
    TEventList *lev = nullptr;

    Int_t kch = 0;
    DENOISE dn;
    WIENER *w;

    Int_t n_points = memorydepth;
    vector<vector<Double_t>> raw;
    vector<vector<Double_t>> wvf;
    vector<TH1D*> haverage;
    vector<TH1D*> hfft;
    Double_t *time = new Double_t[n_points];

    string plot_opt = "AL";
    TGraph *gwvf;
    string xlabel = "Time (ns)";
    string ylabel = "Amplitude (ADC Channels)";

    Double_t ymin = 0;
    Double_t ymax = 0;

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    void setAnalyzer(){
      w = new WIENER(myname.c_str(),dtime,250,1e-9,1e6,memorydepth);
      if(f == nullptr) f = new TFile("analyzed.root","READ");
      if(t1 == nullptr) t1 = (TTree*)f->Get("t1");
      TList *lb = (TList*)t1->GetListOfBranches();
      lev = new TEventList(Form("lev_%s",myname.c_str()),Form("lev_%s",myname.c_str()));
      nchannels = lb->GetEntries();
      b.resize(nchannels);
      schannel.resize(nchannels);
      ch.resize(nchannels);
      raw.resize(nchannels);
      wvf.resize(nchannels);
      haverage.resize(nchannels);
      hfft.resize(nchannels);

      for (Int_t i = 0; i < nchannels; i++) {

        schannel[i] = lb->At(i)->GetName();
        // schannel[i] = Form("Ch%d",channels[i]);
        b[i] = t1->GetBranch(schannel[i].c_str());
        b[i]->SetAddress(&ch[i]);
        raw[i].resize(n_points);
        wvf[i].resize(n_points);
      }
      for (int j = 0; j < n_points; j++) {
        time[j] = j*dtime;
      }
      nentries = t1->GetEntries();
    }

    void getFFT(Double_t *_v = nullptr){
      if(_v == nullptr) _v = ch[kch].wvf;
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,_v[i]);
      }
      w->fft(w->hwvf);
    }

    void setAnalyzerExt(TFile *&ft, TTree *&tr, vector<TBranch *> &bt){
      f = ft;
      tr = (TTree*)ft->Get("t1");
      t1 = tr;
      bt.resize(nchannels);
      setAnalyzer();
      for (Int_t i = 0; i < nchannels; i++) bt[i] = b[i];
    }

    Double_t integrate(Int_t channel = 0, Double_t from = 0, Double_t to = 0){
      Double_t res = 0;
      if (to == 0) to = memorydepth*dtime;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        res += ch[channel].wvf[i];
      }
      return res*dtime;
    }


    void getWaveform(Int_t myevent = 0, Int_t k = 0, Double_t factor = 1){
      if (k>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return;
      }
      b[k]->GetEvent(myevent);
    }

    bool getWaveformHard(Int_t myevent = 0, Double_t factor = 1){
      if (kch>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return false;
      }
      b[kch]->GetEvent(myevent);
      for (int j = 0; j < n_points; j++) {
        raw[kch][j] = ch[kch].wvf[j]*factor;
        wvf[kch][j] = raw[kch][j];
        time[j] = j*dtime;
      }
      return true;

    }


    void applyMovingAverage(Int_t mafilter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      Double_t *_temp = new Double_t[memorydepth];
      if(mafilter!=0) {
        if (_raw == _filtered) {
          _raw = _temp;
          _filtered = ch[kch].wvf;
          for (Int_t i = 0; i < memorydepth; i++) {
            _raw[i] = _filtered[i];
          }
        }
        dn.movingAverage(_raw,_filtered,mafilter);
      }
    }

    void applyDenoise(Int_t filter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      if (_raw == nullptr){
        _raw = ch[kch].wvf;
        _filtered = ch[kch].wvf;
      }
      dn.TV1D_denoise(_raw,_filtered,n_points,filter);
    }

    void setFreqFilter(Double_t frequency_cut, string filter_type = "gaus"){
      w->setFilter(frequency_cut,filter_type);
    }
    void applyFreqFilter(Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch].wvf;
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->apply_filter();
      w->backfft(*w);
      for(Int_t i = 0; i < memorydepth; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
      // w->hwvf->Draw("");
    }

    void applyBandCut(WIENER *_temp = nullptr, Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch].wvf;
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->applyBandCut(_temp);
      w->backfft(*w);
      for(Int_t i = 0; i < memorydepth; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
    }

    void makeCopy(Double_t *cpy, Double_t *original = nullptr){
      if (original == nullptr) original = ch[kch].wvf;
      for (Int_t i = 0; i < memorydepth; i++) {
        cpy[i] = original[i];
      }
    }

    void drawGraph(string opt = "", Int_t n = memorydepth, Double_t* x = nullptr, Double_t* y = nullptr){
      if (opt == "") opt = plot_opt;
      if (x == nullptr) x = time;
      if (y == nullptr) y = ch[kch].wvf;
      gwvf = new TGraph(n,x,y);
      gwvf->Draw(opt.c_str());
      gwvf->GetXaxis()->SetTitle(xlabel.c_str());
      gwvf->GetYaxis()->SetTitle(ylabel.c_str());

      if(ymax!=0 && ymin!=0){
        gwvf->GetYaxis()->SetRangeUser(ymin,ymax);
      }
      gwvf->SetEditable(kFALSE);
    }

    ANALYZER(string m_myname = "z") : myname{m_myname}{

    }

    void averageWaveform(Int_t maxevent = 0, string selection = ""){
      if (maxevent==0) {
        maxevent = nentries;
      }
      haverage[kch] = new TH1D(Form("haverage_%s_Ch%d",myname.c_str(),kch),"Averaged waveform",memorydepth,0,memorydepth*dtime);
      Int_t total = 0;
      t1->Draw(Form(">>lev_%s",myname.c_str()),selection.c_str());
      Int_t nev = lev->GetN();
      if (maxevent < nev) {
        nev = maxevent;
      }
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev);
        total += 1;
        for (Int_t j = 0; j < memorydepth; j++){
          haverage[kch]->AddBinContent(j+1,ch[kch].wvf[j]);
        }

      }

      haverage[kch]->Scale(1./total);
    }




      void showFFT(Int_t naverage, Int_t maxevent, Int_t dt);
      void averageFFT(Int_t maxevent, string selection);


};



void ANALYZER::averageFFT(Int_t maxevent = 0, string selection = ""){
  if (maxevent==0) {
    maxevent = nentries;
  }
  t1->Draw(Form(">>lev_%s",myname.c_str()),selection.c_str());
  Int_t nev = lev->GetN();
  if (maxevent < nev) {
    nev = maxevent;
  }
  Int_t iev = 0;
  hfft[kch] = (TH1D*)w->hfft->Clone("h");
  Int_t total = 0;
  for(Int_t i = 0; i < nev; i++){
    iev = lev->GetEntry(i);
    getWaveform(iev);
    getWaveform(i);
    getFFT();
    for (Int_t j = 0; j < memorydepth/2; j++) hfft[kch]->AddBinContent(j+1,w->hfft->GetBinContent(j+1));
    total++;
  }
  hfft[kch]->Scale(1./total);
}
void ANALYZER::showFFT(Int_t naverage = 10, Int_t maxevent = 0, Int_t dt = 0){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TH1D *h = (TH1D*)w->hfft->Clone("h");
  vector<TH1D *> hfft(naverage);
  Int_t k = 0;
  TCanvas *c1 = new TCanvas("c1");
  for(Int_t i = 0; i < maxevent; i++){
    getWaveform(i);
    getFFT();
    if (i < naverage) {
      hfft[k] = (TH1D*)w->hfft->Clone(Form("h%d",k));
      for (Int_t j = 0; j < memorydepth/2; j++) h->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
      if (naverage == 1) h->Draw();
    }
    else{
      if (k < naverage) {
        for (Int_t j = 0; j < memorydepth/2; j++) {
          h->AddBinContent(j+1,-hfft[k]->GetBinContent(j+1));
          hfft[k]->SetBinContent(j+1,w->hfft->GetBinContent(j+1));
          h->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
        }
        printf("\rEvent %d", i);
        fflush(stdout);
        c1->Modified();
        c1->Update();
        if (gSystem->ProcessEvents())
          break;
        if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
      }
      else{
        if (i == naverage){
          h->Draw();
        }
        k = 0;
      }
    }

    k++;
    if(naverage == 1) k = 0;
  }
}
