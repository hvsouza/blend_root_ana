// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer.C
// Created: 2022
// ________________________________________ //
//
#include "MYCODES.h"



class ANALYZER{
  public:

    TFile *f;
    TTree *t1;
    vector<TBranch*> b;
    vector<ADC_DATA> ch;
    Int_t nchannels = 0;
    vector<Int_t> channels = {1,2};
    vector<string> schannel;
    Double_t dtime = 4;
    string myname;

    DENOISE dn;
    WIENER *w;

    Int_t n_points = memorydepth;
    Double_t **raw;
    Double_t **wvf;
    Double_t *time = new Double_t[n_points];

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    // void setAnalyzer
    void setAnalyzer(){
      w = new WIENER(myname.c_str(),dtime,250,1e-9,1e6,memorydepth);
      if(f == nullptr) f = new TFile("analyzed.root","READ");
      if(t1 == nullptr) t1 = (TTree*)f->Get("t1");
      TList *lb = (TList*)t1->GetListOfBranches();
      nchannels = lb->GetEntries();
      b.resize(nchannels);
      schannel.resize(nchannels);
      ch.resize(nchannels);

      raw = new Double_t *[nchannels];
      wvf = new Double_t *[nchannels];
      for (Int_t i = 0; i < nchannels; i++) {
        schannel[i] = lb->At(i)->GetName();
        // schannel[i] = Form("Ch%d",channels[i]);
        b[i] = t1->GetBranch(schannel[i].c_str());
        b[i]->SetAddress(&ch[i]);
        raw[i] = new Double_t [n_points];
        wvf[i] = new Double_t [n_points];
      }
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
      b[k]->GetEvent(myevent);
      for (int j = 0; j < n_points; j++) {
        raw[k][j] = ch[k].wvf[j]*factor;
        wvf[j][k] = raw[j][k];
        time[j] = j*dtime;
      }
    }


    ANALYZER(string m_myname = "z") : myname{m_myname}{
    }






};
