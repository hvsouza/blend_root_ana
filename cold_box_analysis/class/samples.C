// ________________________________________ //
// Author: Henrique Souza
// Filename: samples.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"



class SAMPLE: public ANALYZER{
  public:
    Int_t n_points = memorydepth;
    string file = "analyzed.root";
    Int_t kch = 0;
    string plot_opt = "ALP";
    string tree = "t1";
    TGraph *gwvf;
    TH2D *hpers;
    string xlabel = "Time (ns)";
    string ylabel = "Amplitude (ADC Channels)";

    Double_t ymin = 0;
    Double_t ymax = 0;



    void print(){
      t1->Print();
    }

    void setChannel(string mych = "Ch1"){
      for(Int_t i = 0; i < nchannels; i++){
        if(mych == schannel[i])
        {
          kch = i;
          return;
        }
      }
    }

    bool getWaveform(Int_t myevent = 0, Double_t factor = 1){
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
    void applyMovingAverage(Int_t mafilter = 0){
      if(mafilter!=0) {
        vector<Double_t> mawvf = dn.movingAverage(&wvf[kch][0],mafilter);
        for (int j = 0; j < n_points; j++) {
          wvf[kch][j] = mawvf[j];
        }
      }
    }

    void applyDenoise(Int_t filter = 0){
      dn.TV1D_denoise(&wvf[kch][0],&wvf[kch][0],n_points,filter);
    }

    void drawGraph(string opt = "", Int_t n = memorydepth, Double_t* x = nullptr, Double_t* y = nullptr){
      if (opt == "") opt = plot_opt;
      if (x == nullptr) x = time;
      if (y == nullptr) y = wvf[kch];
      gwvf = new TGraph(n,x,y);
      gwvf->Draw(opt.c_str());
      gwvf->GetXaxis()->SetTitle(xlabel.c_str());
      gwvf->GetYaxis()->SetTitle(ylabel.c_str());

      if(ymax!=0 && ymin!=0){
        gwvf->GetYaxis()->SetRangeUser(ymin,ymax);
      }
      gwvf->SetEditable(kFALSE);
    }

    void applyFreqFilter(Double_t frequency_cut, string filter_type = "gaus"){
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,wvf[kch][i]);
      }
      w->fft(w->hwvf);
      w->setFilter(frequency_cut,filter_type);
      w->apply_filter();
      w->backfft(*w);
      for(Int_t i = 0; i < memorydepth; i++){
        wvf[kch][i] = w->hwvf->GetBinContent(i+1);
      }
      // w->hwvf->Draw("");
    }
    void sample_plot(Int_t myevent = 0, Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
      bool state = getWaveform(myevent);
      if (!state) return;
      applyMovingAverage(mafilter);
      applyDenoise(filter);

      drawGraph(plot_opt,n_points,&time[0],&wvf[kch][0]);
    }

    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="0==0"){

      hpers = new TH2D("hpers","hpers",n_points,0,n_points*dtime,nbins,ymin,ymax);
      TCanvas *c1 = new TCanvas();
      for (int i=0; i < t1->GetEntries(); i++) {
        if(cut != "0==0"){
          Int_t nselec = t1->Draw(Form("%s.wvf[%d][]:Iteration$*%f",schannel[kch].c_str(),kch,dtime),Form("Entry$==%d && %s",i,cut.c_str()),"goff");
          if(nselec == 0) continue;
          t1->GetSelectedRows();
          raw[kch] = t1->GetV1();
          time = t1->GetV2();
        }
        else{
          b[kch]->GetEvent(i);

          for (int j = 0; j < n_points; j++) {
            raw[kch][j] = ch[kch].wvf[j];
            wvf[kch][j] = raw[kch][j];
            time[j] = j*dtime;
          }

        }
        applyDenoise(filter);

        for (int j = 0; j < n_points; j++) {
          if(filter==0) wvf[kch][j] = raw[kch][j];
          hpers->Fill(time[j],wvf[kch][j]);
        }

      }


      hpers->Draw("colz");
      hpers->GetXaxis()->SetTitle("Time (ns)");
      hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
      
    }
    SAMPLE(string m_myname = "s"){
      myname = m_myname;
      setAnalyzer();
    }
};
