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
    string tree = "t1";
    TH2D *hpers;



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

      printf("%s not found, run s.print() to check the branches\n",mych.c_str());
    }

    bool getWaveform(Int_t myevent = 0, Double_t factor = 1){
      bool status = getWaveformHard(myevent,factor);
      return status;

    }


    void sample_plot(Int_t myevent = 0, string opt = "", Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
      if (opt == "") opt = plot_opt;
      bool state = getWaveform(myevent);
      if (!state) return;
      applyMovingAverage(mafilter);
      applyDenoise(filter);

      drawGraph(opt,n_points,&time[0],&ch[kch].wvf[0]);
    }

    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="0==0"){

      hpers = new TH2D("hpers","hpers",n_points,0,n_points*dtime,nbins,ymin,ymax);
      TCanvas *c1 = new TCanvas();
      
      for (int i=0; i < t1->GetEntries(); i++) {
        if(cut != "0==0"){
          Int_t nselec = t1->Draw("1",Form("Entry$==%d && %s",i,cut.c_str()),"goff");
          if(nselec == 0) continue;
          // t1->GetSelectedRows();
          // Double_t *traw= t1->GetV1();
          // for (int j = 0; j < n_points; j++) {
            // ch[kch].wvf[j] = traw[j];
          // }
          b[kch]->GetEvent(i);
        }
        else{
          b[kch]->GetEvent(i);
        }
        applyDenoise(filter);

        for (int j = 0; j < n_points; j++) {
          hpers->Fill(time[j],ch[kch].wvf[j]);
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
