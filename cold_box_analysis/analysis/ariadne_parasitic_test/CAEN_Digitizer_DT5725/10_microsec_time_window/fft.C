#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

 #include "TH1D.h"
 #include "TVirtualFFT.h"
 #include "TF1.h"
 #include "TCanvas.h"
 #include "TMath.h"

void shift_waveform(vector<Double_t> &y_raw){
  vector<Double_t> ytemp = y_raw; 
  Int_t old_peak_start = 9900/4; // the deconvolved waveform peak starts more or less here.
    Int_t new_peak_start = 4000/4; // I will make it start at 800 ns
  Int_t old_wave_ref = old_peak_start - new_peak_start; // this will be my reference, see the first for
  // starting from old_wave_ref
  for(Int_t i = 0; i<(memorydepth-old_wave_ref); i++){ // so ytemp is not saturated 
    y_raw[i] = ytemp[old_wave_ref+i];
  }
  Int_t aux = 0;
  for(Int_t i = memorydepth-old_wave_ref; i<memorydepth; i++){
    y_raw[i] = ytemp[aux];
    aux++;
  }
    
}

 
 void fft()
 {
   
   
   Int_t channel = 1;
   
   TFile *fwvf = new TFile("led_response.root","READ");
   TFile *fsphe = new TFile("few_pe_response.root","READ");

   TH1D *hsignal = (TH1D*)fwvf->Get("averaged_led");
   TH1D *hsphe = (TH1D*)fsphe->Get("averaged_few_pe");
   
   
   vector<Double_t> x(memorydepth);
   vector<Double_t> yraw(memorydepth);
   vector<Double_t> ysphe(memorydepth);
   vector<Double_t> ysphe2(memorydepth);
   vector<Double_t> ydeconv(memorydepth);

   for(Int_t i = 0; i<memorydepth; i++){
     yraw[i] = hsignal->GetBinContent(i+1);
     ydeconv[i] = yraw[i];
     ysphe2[i] = hsphe->GetBinContent(i+1);
     x[i] = i*4;
   }
   TGraph *graw = new TGraph(memorydepth,&x[0],&yraw[0]);
   TGraph *gsphe = new TGraph(memorydepth,&x[0],&ysphe2[0]);

   
   
 //   TF1 *func_fit2 = new TF1("func_fit2","([0]*exp(-(x-[2])/[1])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([1]*[1])))*TMath::Erfc((([2]-x)/[3]+[3]/[1])/TMath::Power(2,0.5))+([4]*exp(-(x-[2])/[5])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([5]*[5])))*TMath::Erfc((([2]-x)/[3]+[3]/[5])/TMath::Power(2,0.5))+([6]*exp(-(x-[2])/[7])/TMath::Power(2*TMath::Pi(),0.5)*exp(-[3]*[3]/([7]*[7])))*TMath::Erfc((([2]-x)/[3]+[3]/[7])/TMath::Power(2,0.5)) + [8]*TMath::Erfc((([2]-x)/[3])/TMath::Power(2,0.5))",0,20000);
//    func_fit2->SetNpx(20000);
//    // func_fit2->SetParameters(17.43565,245.695,400.557,29.4995,3.93622,245.715,6.72069,159.092,-7);
//    func_fit2->SetParameters(176,909,324,45.7,-156,1027,25.77,88.9,0);

//    // func_fit2->SetParLimits(8,-8,-1);
//    func_fit2->FixParameter(8,0);
//    gsphe->Fit("func_fit2","R0");
// //    TF1 *flar = new TF1("flar","(([0]/[1])*exp(-(x-6000)/[1])+([2]/[3])*exp(-(x-6000)/[3])) * ((x > 6000) ? 1/([0]/[1]+[2]/[3]) : 0);",0,20000);
  
   // Double_t original_offset = func_fit2->GetParameter(2);
   // Double_t offset = 0;
   // func_fit2->SetParameter(2,offset+original_offset);
//    func_fit2->SetParameter(8,0);
   // TF1 *fpol1 = new TF1("fpol1","pol8",0,280);
   // TF1 *fpol2 = new TF1("fpol2","pol8",280,1030);
   // fpol1->SetNpx(10000);
   // fpol2->SetNpx(10000);
   // fpol1->SetParameters(1,1,1,1,1,1,1,1);
   // fpol2->SetParameters(1,1,1,1,1,1,1,1);
   
   // gsphe->Fit("fpol1","R0");
   // gsphe->Fit("fpol2","R0");


   TF1 *flar = new TF1("flar","[0]*exp(-(x-1000)/[1])+[2]*exp(-(x-1000)/[3])",800,memorydepth);
   flar->SetParameters(0.77,20,0.23,1468.59);
   flar->SetNpx(10000);
   
   Double_t factor = 0;
   Double_t smoothness = 10;
   for(Int_t i = 0; i<memorydepth; i++){
     ysphe[i] = hsphe->GetBinContent(i+1);
    
     //or this
     // ysphe[i] = func_fit2->Eval(i*4);
     // hsphe->SetBinContent(i+1,ysphe[i]);

   }
   
   TCanvas *c1 = new TCanvas();
   c1->cd();
   // graw->Draw("ALP");
   
   Int_t n = memorydepth;
 
   hsignal->Draw();
   hsignal->GetXaxis()->SetLabelSize(0.05);
   hsignal->GetYaxis()->SetLabelSize(0.05);

   TCanvas *c2 = new TCanvas();
   c2->cd();

   //Compute the transform and look at the magnitude of the output
   TH1 *hm=0;
   TVirtualFFT::SetTransform(0);
   hm = (TH1*)hsignal->FFT(hm, "MAG R2C measure");
   hm->SetName("fft_signal");
   hm->SetTitle("Magnitude of the 1st transform");
   hm->Draw();
   //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
   //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
   
   hm->SetStats(kFALSE);
   hm->GetXaxis()->SetLabelSize(0.05);
   hm->GetYaxis()->SetLabelSize(0.05);
   
   TCanvas *c3 = new TCanvas();
   c3->cd();
   TGraph *gsphe2 = new TGraph(memorydepth,&x[0],&ysphe2[0]);
   gsphe2->Draw("ALP");
   TGraph *gsphe3 = new TGraph(memorydepth,&x[0],&ysphe[0]);
   gsphe3->SetLineColor(kRed);
   gsphe3->SetLineWidth(2);
   gsphe3->Draw("SAME");


   
   //Look at the DC component and the Nyquist harmonic:
   //That's the way to get the current transform object:
   TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();


   //Use the following method to get the full output:
   Double_t *re_signal = new Double_t[n];
   Double_t *im_signal = new Double_t[n];
   fft->GetPointsComplex(re_signal,im_signal);
   


   TH1 *hm2 =0;
   TVirtualFFT::SetTransform(0);
   hm2 = hsphe->FFT(hm2, "MAG R2C measure");
   hm2->SetTitle("Magnitude of the 1st transform");
   // hm2->Draw();
   //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
   //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!
   
   
   //Look at the DC component and the Nyquist harmonic:
   //That's the way to get the current transform object:
   TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();


   //Use the following method to get the full output:
   Double_t *re_sphe = new Double_t[n];
   Double_t *im_sphe = new Double_t[n];
   fft2->GetPointsComplex(re_sphe,im_sphe);
   
   
   TCanvas *c4 = new TCanvas();
   c4->cd();
   Double_t *re_final = new Double_t[n];
   Double_t *im_final = new Double_t[n];

   Double_t *re_weiner = new Double_t[n];
   Double_t *im_weiner = new Double_t[n];
  
   TF1 *filter = new TF1("filter","TMath::Gaus(x,[0],[1])",0,n);	// A gaussian filter
   Double_t cutoff_frequency = 90; // In MHz
   Double_t adc_frequency = 250; // In MHz
   Double_t converted_freq = cutoff_frequency*memorydepth/adc_frequency/4;//
   filter->SetParameters(0,converted_freq); // USE THIS!!! -> center = 0, cut = 180

    // TF1 *filter = new TF1("filter","exp(-pow(x/[1],[0]))",0,n);	// Wiener-inspered filter
   // filter->SetParameters(3,70); // USE THIS!!! -> center = 0, cut = 180

   
   for(Int_t i = 0; i<n; i++){
     // // this is to remove any offset that remains, but I am not sure if I should use it now
     // if(i==0){
     //   re_final[i] = 0;
     //   im_final[i] = 0;
     //   continue;
     // }
     re_final[i] = filter->Eval(i)*(re_signal[i]*re_sphe[i] + im_signal[i]*im_sphe[i])/(pow(re_sphe[i],2)+pow(im_sphe[i],2));
     im_final[i] = filter->Eval(i)*(re_sphe[i]*im_signal[i] - re_signal[i]*im_sphe[i])/(pow(re_sphe[i],2)+pow(im_sphe[i],2));

     // re_final[i] = re_signal[i];
     // im_final[i] = im_signal[i];


   }
   
     //Now let's make a backward transform:
   TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &n, "C2R M K");
   fft_final->SetPointsComplex(re_final,im_final);
   fft_final->Transform();
   TH1 *hfinal = 0;
   //Let's look at the output
   hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
   hfinal->Scale(1./memorydepth);
   for(Int_t i = 0; i<n; i++){
     ydeconv[i] = hfinal->GetBinContent(i+1);
   }
   shift_waveform(ydeconv);
   // renorm(ydeconv);
   // normalize(ydeconv);
   TGraph *gconv = new TGraph(memorydepth,&x[0],&ydeconv[0]);
   gconv->GetXaxis()->SetTitle("Time (ns)");
   gconv->GetYaxis()->SetTitle("Amplitude (A. U.)");
   gconv->GetYaxis()->SetTitleSize(0.04);
   gconv->GetYaxis()->SetTitleOffset(0.9);
   gconv->Draw("AL");
//    hfinal->SetTitle("The backward transform result");
//    hfinal->Draw("hist");
//    //NOTE: here you get at the x-axes number of bins and not real values
//    //(in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
//    //also here the y-axes has to be rescaled (factor 1/bins)
// //    hfinal->SetStats(kFALSE);
//    hfinal->GetXaxis()->SetLabelSize(0.05);
//    hfinal->GetYaxis()->SetLabelSize(0.05);
   delete fft_final;
   fft_final=0;


   //  //Now let's make a backward transform:
   // TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
   // fft_back->SetPointsComplex(re_signal,im_signal);
   // fft_back->Transform();
   // TH1 *hb = 0;
   // //Let's look at the output
   // hb = TH1::TransformHisto(fft_back,hb,"Re");
   // hb->SetTitle("The backward transform result");
   // // hb->Draw();
   // //NOTE: here you get at the x-axes number of bins and not real values
   // //(in this case 25 bins has to be rescaled to a range between 0 and 4*Pi;
   // //also here the y-axes has to be rescaled (factor 1/bins)
   // hb->SetStats(kFALSE);
   // hb->GetXaxis()->SetLabelSize(0.05);
   // hb->GetYaxis()->SetLabelSize(0.05);
   // delete fft_back;
   // fft_back=0;
   
   delete [] re_signal;
   delete [] im_signal;
   delete [] re_sphe;
   delete [] im_sphe;
 }
 
