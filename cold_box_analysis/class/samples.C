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
    string tree = "t1";



    void print(){
      t1->Print();
    }


    void sample_plot(Int_t myevent = 0, string opt = "", Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
      if (opt == "") opt = plot_opt;
      bool state = getWaveformHard(myevent,factor);
      if (!state) return;
      applyMovingAverage(mafilter);
      applyDenoise(filter);

      drawGraph(opt,n_points,&time[0],&ch[kch].wvf[0]);
    }


    void showWaveform(Int_t maxevent = 0, Int_t filter = 0, Int_t dt = 0){

      if (maxevent==0) {
        maxevent = nentries;
      }
      TCanvas *c1 = new TCanvas("c1");

      for(Int_t i = 0; i < maxevent; i++){
        sample_plot(i,"AL",filter);
        printf("\rEvent %d", i);
        fflush(stdout);
        c1->Modified();
        c1->Update();
        if (gSystem->ProcessEvents())
          break;
        if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
      }
    }

    SAMPLE(string m_myname = "s", string m_filename = "analyzed.root"){
      myname = m_myname;
      setAnalyzer(m_filename);
    }
};
