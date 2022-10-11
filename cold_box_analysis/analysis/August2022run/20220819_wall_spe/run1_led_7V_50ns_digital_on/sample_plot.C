#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"


void sample_plot(Int_t myevent = 0){
 DENOISE dn;

 TFile *f = new TFile("analyzed.root","READ");
 TTree *t1 = (TTree*)f->Get("t1");
 ADC_DATA ch1;
 TBranch *bch = t1->GetBranch("Ch1");
 bch->SetAddress(&ch1);

 vector<Double_t> val(memorydepth);
 vector<Double_t> wvf(memorydepth);
 vector<Double_t> time(memorydepth);

 for (int i=myevent; i < myevent+1; i++) {
   bch->GetEvent(i);

   for (int j = 0; j < memorydepth; j++) {
     val[j] = ch1.wvf[j];
     time[j] = j*4;
   }
 }


 dn.TV1D_denoise(&val[0],&wvf[0],memorydepth,16);

 TGraph *gwvf = new TGraph(memorydepth,&time[0],&wvf[0]);
 gwvf->Draw("ALP");

}
