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

    DENOISE dn;
    WIENER *w;

    Int_t n_points = memorydepth;
    Double_t *raw = new Double_t[n_points];
    Double_t *wvf = new Double_t[n_points];
    Double_t *time = new Double_t[n_points];

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    void setAnalyzer(TFile *&ft, TTree *&tr, vector<TBranch *> &bt){
      if(ft==nullptr){
        ft = new TFile("analyzed.root","READ");
      }
      f = ft;
      tr = (TTree*)ft->Get("t1");
      t1 = tr;

      nchannels = channels.size();
      b.resize(nchannels);
      bt.resize(nchannels);
      schannel.resize(nchannels);
      ch.resize(nchannels);
      for (Int_t i = 0; i < nchannels; i++) {
        schannel[i] = Form("Ch%d",channels[i]);
        bt[i] = t1->GetBranch(schannel[i].c_str());
        bt[i]->SetAddress(&ch[i]);
        b[i] = bt[i];
      }

    }

    Double_t integrate(Int_t channel = 0, Double_t from = 0, Double_t to = 0){
      Double_t res = 0;
      if (to == 0) to = memorydepth*dtime;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        res += ch[channel].wvf[i];
      }
      return res*dtime;
    }

    


    ANALYZER(string myname = "z"){
      w = new WIENER(myname.c_str(),dtime,250,1e-9,1e6,memorydepth);
    }





};
