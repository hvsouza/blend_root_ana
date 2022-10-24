// ________________________________________ //
// Author: Henrique Souza
// Filename: samples.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"



class SAMPLE{
  public:
    Int_t n_points = memorydepth;
    string file = "analyzed.root";
    TFile *f;
    TTree *t1;
    Int_t channel = 1;
    string plot_opt = "ALP";
    string tree = "t1";
    Double_t dtime = 4;
    TGraph *gwvf;
    TH2D *hpers;
    string xlabel = "Time (ns)";
    string ylabel = "Amplitude (ADC Channels)";

    Double_t ymin = 0;
    Double_t ymax = 0;

    DENOISE dn;
    SPHE spe;
    ADC_DATA ch;
    WIENER *w = new WIENER("w",dtime,250,1e-9,1e6,memorydepth);

    Double_t *raw = new Double_t[n_points];
    Double_t *wvf = new Double_t[n_points];
    Double_t *time = new Double_t[n_points];

    void print(){
      f = new TFile(file.c_str(),"READ");
      t1 = (TTree*)f->Get(tree.c_str());
      t1->Print();
    }

    bool getWaveform(Int_t myevent = 0, Double_t factor = 1){
      f = new TFile(file.c_str(),"READ");
      t1 = (TTree*)f->Get(tree.c_str());
      string schannel = Form("Ch%d",channel);
      TBranch *bch = t1->GetBranch(schannel.c_str());
      if (bch == nullptr){
        cout << schannel << " not present in the TTree, execute print() to check channels" << endl;
        return false;
      }
      bch->SetAddress(&ch);

      bch->GetEvent(myevent);
      for (int j = 0; j < n_points; j++) {
        raw[j] = ch.wvf[j]*factor;
        wvf[j] = raw[j];
        time[j] = j*dtime;
      }
      return true;

    }
    void applyMovingAverage(Int_t mafilter = 0){
      if(mafilter!=0) {
        vector<Double_t> mawvf = spe.movingAverage(&wvf[0],mafilter);
        for (int j = 0; j < n_points; j++) {
          wvf[j] = mawvf[j];
        }
      }
    }

    void applyDenoise(Int_t filter = 0){
      dn.TV1D_denoise(&wvf[0],&wvf[0],n_points,filter);
    }

    void drawGraph(string opt = "", Int_t n = memorydepth, Double_t* x = nullptr, Double_t* y = nullptr){
      if (opt == "") opt = plot_opt;
      if (x == nullptr) x = time;
      if (y == nullptr) y = wvf;
      gwvf = new TGraph(n,x,y);
      gwvf->Draw(opt.c_str());
      gwvf->GetXaxis()->SetTitle(xlabel.c_str());
      gwvf->GetYaxis()->SetTitle(ylabel.c_str());

      if(ymax!=0 && ymin!=0){
        gwvf->GetYaxis()->SetRangeUser(ymin,ymax);
      }
      gwvf->SetEditable(kFALSE);
    }

    void applyGausLowPass(Int_t frequency_cut){
      for(Int_t i = 0; i < memorydepth; i++){
        w->hwvf->SetBinContent(i+1,wvf[i]);
      }
      w->fft(w->hwvf);
      w->apply_filter(frequency_cut);
      w->backfft(*w);
      for(Int_t i = 0; i < memorydepth; i++){
        wvf[i] = w->hwvf->GetBinContent(i+1);
      }
      // w->hwvf->Draw("");
    }
    void sample_plot(Int_t myevent = 0, Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
      bool state = getWaveform(myevent);
      if (!state) return;
      applyMovingAverage(mafilter);
      applyDenoise(filter);

      drawGraph(plot_opt,n_points,&time[0],&wvf[0]);
    }

    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="0==0"){

      f = new TFile(file.c_str(),"READ");
      t1 = (TTree*)f->Get(tree.c_str());
      ADC_DATA ch;
      string schannel = Form("Ch%d",channel);
      TBranch *bch = t1->GetBranch(schannel.c_str());
      bch->SetAddress(&ch);

      hpers = new TH2D("hpers","hpers",n_points,0,n_points*dtime,nbins,ymin,ymax);
      TCanvas *c1 = new TCanvas();
      for (int i=0; i < t1->GetEntries(); i++) {
        if(cut != "0==0"){
          Int_t nselec = t1->Draw(Form("%s.wvf[]:Iteration$*%f",schannel.c_str(),dtime),Form("Entry$==%d && %s",i,cut.c_str()),"goff");
          if(nselec == 0) continue;
          t1->GetSelectedRows();
          raw = t1->GetV1();
          time = t1->GetV2();
        }
        else{
          bch->GetEvent(i);

          for (int j = 0; j < n_points; j++) {
            raw[j] = ch.wvf[j];
            wvf[j] = raw[j];
            time[j] = j*dtime;
          }

        }
        applyDenoise(filter);

        for (int j = 0; j < n_points; j++) {
          if(filter==0) wvf[j] = raw[j];
          hpers->Fill(time[j],wvf[j]);
        }

      }


      hpers->Draw("colz");
      hpers->GetXaxis()->SetTitle("Time (ns)");
      hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

    }
};
