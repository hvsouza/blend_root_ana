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
    string filename = "analyzed.root";
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
    TH2D *hpers;

    Double_t ymin = 0;
    Double_t ymax = 0;
    Double_t temp_charge = 0;
    Double_t temp_max = 0;

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    void setAnalyzer(string m_filename = "analyzed.root"){
      filename = m_filename;
      w = new WIENER(myname.c_str(),dtime,250,1e-9,1e6,memorydepth);
      if(f == nullptr) f = new TFile(filename.c_str(),"READ");
      if(t1 == nullptr) t1 = (TTree*)f->Get("t1");
      TList *lb = (TList*)t1->GetListOfBranches();
      lev = new TEventList(Form("lev_%s",myname.c_str()),Form("lev_%s",myname.c_str()));
      nchannels = lb->GetEntries();
      b.resize(nchannels);
      schannel.resize(nchannels);
      channels.resize(nchannels);
      ch.resize(nchannels);
      raw.resize(nchannels);
      wvf.resize(nchannels);
      haverage.resize(nchannels);
      hfft.resize(nchannels);

      for (Int_t i = 0; i < nchannels; i++) {

        schannel[i] = lb->At(i)->GetName();
        channels[i] = schannel[i][2] - '0';
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

    void setAnalyzerExt(TFile *&ft, TTree *&tr, vector<TBranch *> &bt){
      f = ft;
      tr = (TTree*)ft->Get("t1");
      t1 = tr;
      bt.resize(nchannels);
      setAnalyzer();
      for (Int_t i = 0; i < nchannels; i++) bt[i] = b[i];
    }

    void getFFT(Double_t *_v = nullptr){
      if(_v == nullptr) _v = ch[kch].wvf;
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,_v[i]);
      }
      w->fft(w->hwvf);
    }

    Double_t getMean(Double_t from, Double_t to){
      Double_t res = 0;
      Double_t totalev = 0;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        res += ch[kch].wvf[i];
        totalev+=1;
      }
      return res/totalev;
    }

    Double_t getMaximum(Double_t from, Double_t to){
      Double_t max = -1e12;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        if(ch[kch].wvf[i]>=max){
          max = ch[kch].wvf[i];
        }
      }
      return max;
    }

    Double_t rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}){
      kch = channel;
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);
      Double_t peak_level = 0;
      if (ispulse){
        peak_level = getMaximum(peak_range[0],peak_range[1]);
      }
      else{
        peak_level = getMean(peak_range[0],peak_range[1]);
      }
      Int_t startpoint = baseline_range[1]/dtime;
      Int_t endpoint = peak_range[1]/dtime;
      Bool_t triggered = false;
      Double_t time_mark = 0;

      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch].wvf[i] >= 0.1*(peak_level-baseline_level) && triggered==false) {
          triggered = true;
          time_mark = i*dtime;
        }
        else if(triggered == true && ch[kch].wvf[i]>=0.9*(peak_level - baseline_level)){
          time_mark = i*dtime - time_mark;
          return time_mark;
        }
      }
      return 0;
    }

    void gen_rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, Double_t filter = 0, TH1D *htemp = nullptr){
      kch = channel;
      for(Int_t i = 0; i < nentries; i++){
        getWaveform(i,kch);
        applyDenoise(filter);
        htemp->Fill(rise_time(kch, baseline_range, ispulse, peak_range));
      }
    }


    Double_t getTOT(Int_t channel = 0, Double_t from = 0, Double_t to = 0, Double_t time_trigger = 20, Double_t minimum_gap_ns = 12){
      Bool_t triggered = false;
      Double_t time_mark = 0;
      for(Int_t i = from/dtime; i<to/dtime; i++){
        if(ch[channel].wvf[i]>time_trigger && triggered == false){
          triggered = true;
          time_mark = i*dtime;
          // break;
        }
        else if(triggered==true && ch[channel].wvf[i]<time_trigger && (i*dtime-time_mark)>=minimum_gap_ns){
          time_mark = i*dtime - time_mark;
          return time_mark;
        }
      }
      return 0;
    }


    void integrate(Int_t channel = 0, Double_t from = 0, Double_t to = 0){
      Double_t res = 0;
      if (to == 0) to = memorydepth*dtime;
      Double_t max = -1e12;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        res += ch[channel].wvf[i];
        if(ch[channel].wvf[i]>=max){
          max = ch[channel].wvf[i];
        }
      }
      temp_charge = res*dtime;
      temp_max = max;
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
      if (filter == 0) return;
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

    void getBackFFT(Double_t *_filtered = nullptr){
      if (_filtered == nullptr){
        _filtered = ch[kch].wvf;
      }
      w->backfft(*w);
      for(Int_t i = 0; i < memorydepth; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }

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


    void selectByAmplitude(Double_t filter = 0, Double_t xmin = 0, Double_t xmax = 0, Double_t limit = 100){
      lev = (TEventList*)gDirectory->Get(Form("lev_%s",myname.c_str()));
      if(lev->GetN() == 0){
        getSelection(myname, "");
      }
      TEventList *ttemp = new TEventList("ttemp", "ttemp");

      if (xmax == 0) {
        xmax = memorydepth*dtime;
      }
      for(Int_t i = 0; i < lev->GetN(); i++){
        getWaveform(lev->GetEntry(i), kch);
        applyDenoise(filter);
        for(Int_t j = xmin/dtime; j < xmax/dtime; j++){
          if(ch[kch].wvf[j] > limit){
            ttemp->Enter(i);
            break;
          }
        }
      }
      lev->Subtract(ttemp);
      delete ttemp;
    }
    void getSelection(string myname, string selection)
    {
      if(selection!="use_mine"){
        t1->Draw(Form(">>lev_%s",myname.c_str()),selection.c_str());
        lev = (TEventList*)gDirectory->Get(Form("lev_%s",myname.c_str()));
      }
      else{
        if (lev->GetN() == 0){
          cout << "No event selectected... " << endl;
        }
      }
    }

    void averageWaveform(Int_t maxevent = 0, string selection = "", Double_t filter =  0){
      if (maxevent==0) {
        maxevent = nentries;
      }
      haverage[kch] = new TH1D(Form("haverage_%s_Ch%d",myname.c_str(),kch),"Averaged waveform",memorydepth,0,memorydepth*dtime);
      Int_t total = 0;
      getSelection(myname, selection);
      Int_t nev = lev->GetN();
      if (maxevent < nev) {
        nev = maxevent;
      }
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev);
        applyDenoise(filter);
        total += 1;
        for (Int_t j = 0; j < memorydepth; j++){
          haverage[kch]->AddBinContent(j+1,ch[kch].wvf[j]);
        }

      }

      haverage[kch]->Scale(1./total);
    }

    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="", Double_t factor = 1){

      hpers = new TH2D("hpers","hpers",n_points,0,n_points*dtime,nbins,ymin,ymax);
      TCanvas *c1 = new TCanvas();

      getSelection(myname, cut);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        b[kch]->GetEvent(iev);
        applyDenoise(filter);

        for (int j = 0; j < n_points; j++) {
          hpers->Fill(time[j],ch[kch].wvf[j]*factor);
        }

      }


      hpers->Draw("colz");
      hpers->GetXaxis()->SetTitle("Time (ns)");
      hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

    }

    void setChannel(string mych = "Ch1"){
      for(Int_t i = 0; i < nchannels; i++){
        if(mych == schannel[i])
        {
          kch = i;
          return;
        }
      }

      printf("%s not found, run s.print() to check the branches\n",mych.c_str());
    }


    void showFFT(Int_t naverage, Int_t maxevent, Int_t dt, bool inDecibel);
    void averageFFT(Int_t maxevent, string selection, bool inDecibel, Double_t filter);


};



void ANALYZER::averageFFT(Int_t maxevent = 0, string selection = "", bool inDecibel = false, Double_t filter = 0){
  if (maxevent==0) {
    maxevent = nentries;
  }
  getSelection(myname, selection);
  Int_t nev = lev->GetN();
  if (maxevent < nev) {
    nev = maxevent;
  }
  Int_t iev = 0;
  hfft[kch] = (TH1D*)w->hfft->Clone(Form("hfft_%s_ch%d",myname.c_str(),kch));
  hfft[kch]->Reset();
  Int_t total = 0;
  for(Int_t i = 0; i < nev; i++){
    iev = lev->GetEntry(i);
    getWaveform(iev,kch);
    applyDenoise(filter);
    getFFT();
    for (Int_t j = 0; j < memorydepth/2; j++) hfft[kch]->AddBinContent(j+1,w->hfft->GetBinContent(j+1));
    total++;
  }
  hfft[kch]->Scale(1./total);
  hfft[kch]->SetEntries(total);

  if (inDecibel){
    w->convertDecibel(hfft[kch]);
    hfft[kch]->GetYaxis()->SetTitle("Magnitude (dB)");
  }
}
void ANALYZER::showFFT(Int_t naverage = 10, Int_t maxevent = 0, Int_t dt = 0, bool inDecibel = false){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TH1D *h = (TH1D*)w->hfft->Clone("h");
  vector<TH1D *> hfft(naverage);
  Int_t k = 0;
  TCanvas *c1 = new TCanvas("c1");
  for(Int_t i = 0; i < maxevent; i++){
    getWaveform(i,kch);
    getFFT();
    if (inDecibel) w->convertDecibel();
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
