#define memorydepth 250000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void FFT(){
    vector<string> runs = {"run0_light_leakage_test_led_flange_covered","run1_light_leakage_test_led_fibers_connected","run2_light_leakage_test_led_fibers_connected_dcem_digital_off"};

    vector<string> desc = {"DCem digital on (led fibers disconnected)","DCem digital on","DCem Digital off"};
    vector<Color_t> mycolors = {kViolet,kRed,kBlack};

     vector<Int_t> channels = {1};
     vector<string> names = {"miniArapuca 37V (A1ch1)"};
     Int_t n = channels.size();
     vector<ADC_DATA> ch(n);
     vector<TBranch*> bch(n);
     vector<TCanvas*> c(n);

    vector<vector<WIENER*>> wres(3,vector<WIENER*>(n,NULL));
     for(Int_t i = 0; i < n; i++){
         c[i] = new TCanvas(names[i].c_str(),names[i].c_str());
     }
     for (int jrun = 0; jrun < 3; jrun++) {

        TFile *f = new TFile(Form("%s/analyzed.root",runs[jrun].c_str()),"READ");
        TTree *t1 = (TTree*)f->Get("t1");

        Double_t dtime = 2;
        Double_t freq = 500;


        // vector<TH1D*> hres(n);
        for(Int_t i = 0; i < n; i++){
            bch[i] = t1->GetBranch(Form("Ch%i",channels[i]));
            bch[i]->SetAddress(&ch[i]);
            // hres[i] = new TH1D(Form("hres%d",channels[i]), Form("hres%d",channels[i]),memorydepth, 0, dtime * memorydepth);
        }

        TH1D *hsignal = new TH1D("hsignal", "hsignal", memorydepth, 0, dtime * memorydepth);
        WIENER ws("ws",dtime,freq,1e-9,1e6,memorydepth);
        for(Int_t i = 0; i < n; i++){
            wres[jrun][i] = new WIENER(Form("run%d_miniArapuca 37V (A1ch%d)",jrun,i),dtime,freq,1e-9,1e6,memorydepth);
            // for(Int_t k = 0; k < 1; k++){
            for(Int_t k = 0; k < t1->GetEntries(); k++){
                bch[i]->GetEvent(k);
                for (Int_t j = 0; j < memorydepth; j++) {
                    hsignal->SetBinContent(j+1,ch[i].wvf[j]);
                }
                ws.fft(hsignal);
                for (Int_t j = 0; j < memorydepth/2; j++) {
                    wres[jrun][i]->hfft->AddBinContent(j+1,ws.hfft->GetBinContent(j+1));
                }
            }
            wres[jrun][i]->hfft->ResetStats();
            wres[jrun][i]->hfft->SetLineColor(mycolors[jrun]);
            wres[jrun][i]->hfft->SetTitle(desc[jrun].c_str());
            wres[jrun][i]->hfft->Draw("SAME");
            c[i]->cd();
        }


    }
}
