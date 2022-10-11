#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void FFT(){

    TFile *f = new TFile("analyzed.root","READ");
    TTree *t1 = (TTree*)f->Get("t1");
    vector<Int_t> channels = {1};
    vector<string> names = {"miniArapuca 37V (A1ch1)"};
    Int_t n = channels.size();
    vector<ADC_DATA> ch(n);
    vector<TBranch*> bch(n);
    vector<TCanvas*> c(n);

    Double_t dtime = 2;
    Double_t freq = 500;


    // vector<TH1D*> hres(n);
    for(Int_t i = 0; i < n; i++){
        bch[i] = t1->GetBranch(Form("Ch%i",channels[i]));
        bch[i]->SetAddress(&ch[i]);
        c[i] = new TCanvas(names[i].c_str(),names[i].c_str());
        // hres[i] = new TH1D(Form("hres%d",channels[i]), Form("hres%d",channels[i]),memorydepth, 0, dtime * memorydepth);
    }

    TH1D *hsignal = new TH1D("hsignal", "hsignal", memorydepth, 0, dtime * memorydepth);
    WIENER ws("ws",dtime,freq,1e-9,1e6,memorydepth);
    vector<WIENER*> wres(n);
    for(Int_t i = 0; i < n; i++){
        wres[i] = new WIENER(names[i].c_str(),dtime,freq,1e-9,1e6,memorydepth);
        // for(Int_t k = 0; k < 1000; k++){
        for(Int_t k = 0; k < t1->GetEntries(); k++){
            bch[i]->GetEvent(k);
            for (Int_t j = 0; j < memorydepth; j++) {
                hsignal->SetBinContent(j+1,ch[i].wvf[j]);
            }
            ws.fft(hsignal);
            for (Int_t j = 0; j < memorydepth/2; j++) {
                wres[i]->hfft->AddBinContent(j+1,ws.hfft->GetBinContent(j+1));
            }
        }
        wres[i]->hfft->ResetStats();
        c[i]->cd();
        wres[i]->hfft->Draw("");
    }




}
