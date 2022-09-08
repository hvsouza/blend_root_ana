#include "MYCODES.h"


class SAMPLE{
  public:
    Int_t n_points = memorydepth;
    string file = "analyzed.root";
    string channel = "Ch1";
    string plot_opt = "ALP";
    string tree = "t1";
    Double_t dtime = 4;
    TGraph *gwvf;
    TH2D *hpers;

    void sample_plot(Int_t myevent = 0, Int_t filter = 0, Double_t factor = 1., Int_t mafilter = 0){
      DENOISE dn;
      SPHE spe;
      TFile *f = new TFile(file.c_str(),"READ");
      TTree *t1 = (TTree*)f->Get(tree.c_str());
      ADC_DATA ch;
      TBranch *bch = t1->GetBranch(channel.c_str());
      bch->SetAddress(&ch);

      vector<Double_t> val(n_points);
      vector<Double_t> wvf(n_points);
      vector<Double_t> time(n_points);

      bch->GetEvent(myevent);

      for (int j = 0; j < n_points; j++) {
        val[j] = ch.wvf[j]*factor;
        wvf[j] = val[j];
        time[j] = j*dtime;
      }


      if(mafilter!=0) {
        vector<Double_t> mawvf = spe.movingAverage(val,mafilter);
        for (int j = 0; j < n_points; j++) {
          wvf[j] = val[j] = mawvf[j];
        }
      }
      if(filter!=0) dn.TV1D_denoise(&val[0],&wvf[0],n_points,filter);

      gwvf = new TGraph(n_points,&time[0],&wvf[0]);
      gwvf->Draw(plot_opt.c_str());
      gwvf->GetXaxis()->SetTitle("Time (ns)");
      gwvf->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

    }
    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Int_t filter = 0, string cut="0==0"){
      DENOISE dn;

      TFile *f = new TFile(file.c_str(),"READ");
      TTree *t1 = (TTree*)f->Get(tree.c_str());
      ADC_DATA ch;
      TBranch *bch = t1->GetBranch(channel.c_str());
      bch->SetAddress(&ch);

      Double_t *val = new Double_t[n_points];
      Double_t *wvf = new Double_t[n_points];
      Double_t *time = new Double_t[n_points];

      hpers = new TH2D("hpers","hpers",n_points,0,n_points*dtime,nbins,ymin,ymax);
      TCanvas *c1 = new TCanvas();
      for (int i=0; i < t1->GetEntries(); i++) {
        if(cut != "0==0"){
          Int_t nselec = t1->Draw(Form("%s.wvf[]:Iteration$*%f",channel.c_str(),dtime),Form("Entry$==%d && %s",i,cut.c_str()),"goff");
          if(nselec == 0) continue;
          t1->GetSelectedRows();
          val = t1->GetV1();
          time = t1->GetV2();
        }
        else{
          bch->GetEvent(i);

          for (int j = 0; j < n_points; j++) {
            val[j] = ch.wvf[j];
            wvf[j] = val[j];
            time[j] = j*dtime;
          }

        }
        if(filter!=0) dn.TV1D_denoise(&val[0],&wvf[0],n_points,filter);

        for (int j = 0; j < n_points; j++) {
          if(filter==0) wvf[j] = val[j];
          hpers->Fill(time[j],wvf[j]);
        }

      }


      hpers->Draw("colz");
      hpers->GetXaxis()->SetTitle("Time (ns)");
      hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

    }
};
