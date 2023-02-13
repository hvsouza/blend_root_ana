// ________________________________________ //
// Author: Henrique Souza
// Filename: calibrationCodes.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"


class  MyFunctionObject{
  public:

    Int_t n_peaks;
    // use constructor to customize your function object
    Double_t operator()(Double_t *x, Double_t *par) {
      Double_t f;
      Double_t xx = x[0];
      f  = abs(par[0])*exp(-0.5*TMath::Power((xx-par[1])/par[2],2)); // first argument
      f = f+abs(par[3])*exp(-0.5*TMath::Power((xx-par[4])/par[5],2));
      f = f+abs(par[6])*exp(-0.5*TMath::Power((xx-par[7])/(TMath::Power((2),0.5)*par[5]),2));
      for(Int_t i = 1; i<n_peaks; i++){
        f = f+ abs(par[i+7])*exp(-0.5*TMath::Power((xx-(par[4]+(i+1)*(par[7]-par[4])))/(TMath::Power((i+2),0.5)*par[5]),2));
      }
      return f;
    }
};


class  MyFunctionFree{
  public:

    Int_t n_peaks;
    // use constructor to customize your function object
    Double_t operator()(Double_t *x, Double_t *par) {
      Double_t f;
      Double_t xx = x[0];
      f  = abs(par[0])*exp(-0.5*TMath::Power((xx-par[1])/par[2],2)); // first argument
      f = f+abs(par[3])*exp(-0.5*TMath::Power((xx-par[4])/par[5],2));
      f = f+abs(par[6])*exp(-0.5*TMath::Power((xx-par[7])/par[8],2));
      for(Int_t i = 1, j = 1; i<n_peaks; i++){
        f = f+ abs(par[j+8])*exp(-0.5*TMath::Power((xx-(par[4]+(i+1)*(par[7]-par[4])))/par[j+8+1],2));
        j+=2;
      }
      return f;
    }
};



































class Calibration

{
    
    
    
  public:
    


    string myname = "c";
    Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits
    Int_t linhasEvento = 9000;
  
    Int_t channel = 1;
    // _______________ Parameters for giveMeSphe_MeanPulseVersion _______________/
    Double_t time_sample = 6000; // in ns
  
    Double_t timeLow = 3650; // integration limit
    Double_t timeHigh = 4200;
  
    Int_t fn = 0;
  
  
    // _______________ Parameters for fit_sphe_wave _______________/
  
    string rootFile = "";
    string histogram_name = "";
  
    Int_t n_peaks = 7;
    Double_t peak0 = 0.1;
    Double_t mean0 = -900;
    Double_t sigma0 = 500;
  
    Double_t peak1 = 0.004;
    Double_t mean1 = 5000;
    Double_t sigma1 = 400;
  
    Double_t startpoint = 0.001;
    Double_t poisson_ratio = 3.;
    TF1 *faux = nullptr;
  
    Double_t xmin = -10000;
    Double_t xmax = 40000;

    Double_t deltaplus=1.4;
    Double_t deltaminus=1.2;

    Double_t stdVar = 0.5;
  
    Int_t rebin = 4;
  
    Bool_t fixZero = false;
  
    Bool_t make_free_stddevs = false;

    Double_t sphe_charge = 0; // wave0
    Double_t sphe_charge2 = 0; // wave0
  
  
    // This was add here to try fitting dark noise data with already found values
    Bool_t darknoise = false;
  
    Bool_t is_poisson_test = false; // if running tests of poisson statistics
  
  
  
  
  
    // _______________ Parameters for smoothing  _______________/
    Int_t smooth_factor = 35;


    // ____________________________________________________________________________________________________ //
    void fit_sphe_wave(string name, bool optimize = true){

      if (optimize == true) searchParameters(name,2,false);
      makeSphe(name);
    }
    Double_t fact(Double_t n){
      if( n == 1 || n == 0 ){
        return 1;
      }
      else{
        return n*fact(n-1);
      }
    }
    Double_t poisson(Double_t lambda, Double_t k){
      Double_t factk = fact(k);
      Double_t nume = exp(-lambda)*TMath::Power(lambda,k);
      return nume/factk;
    }

    void searchParameters(string histogram, Double_t sigmaSearch = 2, bool debug = false){

      TFile *f1 = nullptr;
      TH1D *h = nullptr;
      if (rootFile != "") {
        f1 = new TFile(rootFile.c_str(),"READ");
        h = (TH1D*)f1->Get(histogram.c_str());
      }
      else{
        h = (TH1D*)gDirectory->Get(histogram.c_str());
      }
      histogram_name = histogram;


      h->Rebin(rebin);
      Double_t scale = 1/(h->Integral());
      h->Scale(scale);
      Int_t nbins = h->GetNbinsX();
      Double_t *source = new Double_t[nbins];
      Double_t *destVector = new Double_t[nbins];

      Double_t sigma;
      Int_t fNPeaks = 0;
      Int_t nfound = 0;
      Int_t bin;
      // const Int_t nbins = 1024;
      xmin = h->GetXaxis()->GetXmin();
      xmax = h->GetXaxis()->GetXmax();
      Double_t a;
      h->SetTitle("High resolution peak searching, number of iterations = 3");
      h->GetXaxis()->SetRange(1,nbins);
      TH1D *d = new TH1D("d","",nbins,xmin,xmax);

      for (Int_t i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);


      TSpectrum *s = new TSpectrum();

      Double_t sigmastep = 0.2;
      while(nfound < 3 && sigmaSearch>=1){
        nfound = s->SearchHighRes(source, destVector, nbins, sigmaSearch, 2, kFALSE, 5, kTRUE, 3);
        sigmaSearch-=sigmastep;
      }
      if(nfound < 3){
        cout << "Could not optimize parameters, only " << nfound << " peaks found" << endl;
        debug = true;
      }
      Double_t *xpeaks_t = s->GetPositionX();
      vector<Double_t> xpeaks(nfound);
      Int_t pos0 = 0; // this in case there is some peak before zero
      vector<Double_t> fPositionX(nfound-pos0);
      vector<Double_t> fPositionY(nfound-pos0);
      for (Int_t i = 0; i < nfound; i++) xpeaks[i] = xpeaks_t[i];
      std::sort(xpeaks.begin(),xpeaks.end());
      for (Int_t i = 0; i < nfound-pos0; i++) {
        a=xpeaks[i+pos0];
        bin = 1 + Int_t(a + 0.5);
        fPositionX[i] = h->GetBinCenter(bin);
        fPositionY[i] = h->GetBinContent(bin);
      }


      Double_t upplim_gaus_base = fPositionX[0]+(3./4)*(fPositionX[1]-fPositionX[0])/2;
      faux = new TF1("faux","gaus(0)",xmin,upplim_gaus_base);
      faux->SetParameters(fPositionY[0],fPositionX[0],sqrt(fPositionX[0]));
      h->Fit("faux","R0Q");

      peak0 = fPositionY[0];
      mean0 = fPositionX[0];
      sigma0 = faux->GetParameter(2);

      peak1 = fPositionY[1];
      mean1 = fPositionX[1];
      sigma1 = 0.8*faux->GetParameter(2);

      Int_t bin_second = 0; // corresponding bin of the second peak
      if (nfound-pos0 >= 3) // in case it was found
      {
        startpoint =  fPositionY[2];
        bin_second = 1 + xpeaks[pos0 + 2] + 0.5;
      }
      else{
        Int_t bin_first_peak = 1 + Int_t(xpeaks[1+pos0]);
        Int_t bin_baseline = 1 + Int_t(xpeaks[0+pos0]);
        bin_second = 2*bin_first_peak-bin_baseline; // same as b + (b-a)
        startpoint = h->GetBinContent(bin_second);
        // cout << startpoint << " " << bin_first_peak << " " << bin_baseline << " " << bin_second << " " << h->GetBinCenter(bin_second) << endl;
      }
      Double_t lowestpt = 0;
      for (Int_t i = bin_second; i < nbins; i++){
        if(h->GetBinContent(i+1)<= 10*h->GetMinimum(0)*scale){
          lowestpt = h->GetBinCenter(i);
          break;
        }
        if (i == nbins-1){
          lowestpt = h->GetBinCenter(i);
        }
      }
      n_peaks = Int_t(lowestpt/(fPositionX[1]-fPositionX[0])); // this should probably be -1 !!
      if(n_peaks <= 0 || nfound-pos0 < 3){
        cout << "Something wrong: npeaks = " << n_peaks << endl;
        cout << "Changing to 5" << endl;
        n_peaks = 5;
      }
      else{
        if(n_peaks > 10){ // this value was 18 at some point
          n_peaks = 10;
        }
        xmax = (n_peaks+1)*(fPositionX[1]-fPositionX[0]) + fPositionX[0];
      }


      Double_t total_events = h->Integral("width");
      Double_t zero_events = faux->Integral(xmin, xmax);
      Double_t prob_zeros = zero_events/total_events;

      Double_t width = h->GetBinWidth(1);
      Double_t startbin = h->GetBinCenter(1);

      // bin_center = (nbins-1)*width + startbin;
      Int_t npois = 3;
      Int_t bin_n_peak = round(npois*mean1 - startbin)/width + 1;

      Double_t lambda1 = -TMath::Log(prob_zeros);
      Double_t lambda2 = h->GetMean()/mean1;
      Double_t lambda = lambda1;

      Double_t peak_by_formula1 = startpoint/(pow(poisson(lambda1,npois-1)/poisson(lambda1,npois),npois-2));
      Double_t peak_by_formula2 = startpoint/(pow(poisson(lambda2,npois-1)/poisson(lambda2,npois),npois-2));
      Double_t diff1 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula1);
      Double_t diff2 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula2);


      lambda = diff1 <= diff2 ? lambda1 : lambda2;

      Double_t poisson2 = poisson(lambda,2);
      Double_t poisson3 = poisson(lambda,3);
      Double_t poisson_ratio = poisson2/poisson3;


      if(debug)
      {
        TCanvas *cdb = new TCanvas("cdb");
        cdb->cd();
        h->GetYaxis()->SetTitle("Normalized count");
        h->GetXaxis()->SetTitle("Charge (ADC*ns)");
        h->Draw("hist");
        TPolyMarker * pm = new TPolyMarker(nfound-pos0, &fPositionX[0], &fPositionY[0]);
        pm->SetMarkerStyle(23);
        pm->SetMarkerColor(kRed);
        pm->SetMarkerSize(1.3);
        pm->Draw("SAME");

        for (Int_t i = 0; i < nbins; i++) d->SetBinContent(i + 1,destVector[i]);
        d->SetLineColor(kRed);
        d->Draw("SAME");

        printf("Found %d candidate peaks\n",nfound);
        for(Int_t i=0;i<nfound;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
        faux->Draw("SAME");
        cout << "npeaks = " << n_peaks << " lowest = " << lowestpt << " spe = " << (fPositionX[1]-fPositionX[0]) << endl;
      }



      // delete f1, h, source, destVector;
    }



    // ____________________________________________________________________________________________________ //
    string startingPump(){
      string f = "gaus(0) + gaus(3)";
      Int_t aux = 0;
      for(Int_t i = 0; i<n_peaks; i++){
        f = f + " + gaus(" + to_string(i+6+aux) + ")";
        aux = aux+2;
      }
      return f;
      //         TF1 *func = new TF1("func","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",-2000,5000);

    }

    // ____________________________________________________________________________________________________ //
    void getMyParameters(Double_t peaks[],Double_t stdpeaks[],TF1 *func){
      for(Int_t i=0; i<n_peaks; i++){
        peaks[i] = (i+2)*(func->GetParameter(4));
        stdpeaks[i] = TMath::Power((i+2),0.5)*func->GetParameter(5);
        //         cout << peaks[i] << endl;
      }
    }

    // ____________________________________________________________________________________________________ //
    void makeSphe(string histogram){

      string histogram_tempo = histogram+"_stat";
      TH1D *hcharge;
      TFile *f1 = nullptr;
      if (rootFile != "") {
        TFile *f1 = new TFile(rootFile.c_str(),"READ");
        hcharge = (TH1D*)f1->Get(histogram.c_str());
      }
      else{
        hcharge = (TH1D*)gDirectory->Get(histogram.c_str());
      }
      histogram_name = histogram;
      hcharge->Rebin(rebin);
    
      ofstream out;
      out.open("sphe.txt", ios::app);

      // ____________________________ Start of sphe fit ____________________________ //
      hcharge->GetYaxis()->SetTitle("Normalized count");
      hcharge->GetYaxis()->SetTitleOffset(1.0);
      hcharge->GetXaxis()->SetTitle("Charge (ADC*nsec)");
    
    
      TCanvas *c1 = new TCanvas("c1","Carga");
      // c1->SetLogy();
      gPad->SetGrid(1,1);
      gPad->SetTicks(1,1);
      gStyle->SetOptFit();
    
      //First function, will almost fit freely
      TF1 *func = new TF1("func",startingPump().c_str(),xmin,xmax);
      TF1 *fu[2+n_peaks];
      string funame;
      for(Int_t i = 0; i<(2+n_peaks); i++){
        funame = "fu_"+to_string(i);
        fu[i] = new TF1(funame.c_str(),"gaus(0)",xmin,xmax);
      }
    
      Double_t peaks[n_peaks];
      Double_t stdpeaks[n_peaks];
    
      if(peak0==0){
        fixZero = true;
      }
    
      func->SetParameters(peak0,mean0,sigma0,peak1,mean1,sigma1); // this values can change
      // func->SetParLimits(2,0.9*sigma0,2*sigma0);
    
      Int_t aux = 0;
      // TO BE CHECKED ! Maybe it is better to not control the startpoing so hardly
      Double_t temp_startpoint = startpoint;
      for(Int_t i = 0; i<n_peaks; i++){
        func->SetParameter((i+6+aux),temp_startpoint);
        aux++;
        func->SetParameter((i+6+aux),(i+2)*(mean1 - mean0) + mean1);
        aux++;
        func->SetParameter((i+6+aux),sqrt(i+2)*sigma1);
        temp_startpoint = temp_startpoint/poisson_ratio;
      }
      func->SetParName(4,"#mu");
      func->SetParName(5,"#sigma");
    
      func->SetNpx(1000);
    
      if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
      }
    
      if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
      }
    
      getMyParameters(peaks,stdpeaks,func);
    
    
      Double_t scale = 1/(hcharge->Integral());
      hcharge->Scale(scale);
      hcharge->Draw("hist");
      hcharge->Fit("func","R0Q");
      // Debug level:
      // 0 none
      // 1 first general fit
      // 2 second general fit
      // 3 first lastOne fit wiht fix parameters

      Int_t debug_level = 0;
      if (debug_level == 1){
        func->Draw("SAME");
        return;
      }

      // recovery parameters
      getMyParameters(peaks,stdpeaks,func);
      aux=0;
      // use as fixed
      for(Int_t i = 0; i<n_peaks; i++){

        func->FixParameter((i+7+aux),peaks[i]-(i+1)*func->GetParameter(1));
        aux++;
        func->FixParameter((i+7+aux),stdpeaks[i]);
        aux++;
      }

      if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
      }

      if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
      }
      hcharge->Fit("func","R0Q");
      if (debug_level == 2){
        func->Draw("SAME");
        return;
      }

      // set a new function, now fixed for real

      MyFunctionObject MyFunc;
      MyFunc.n_peaks = n_peaks;

      TF1 *lastOne = new TF1("lastOne",MyFunc,xmin,xmax,7+n_peaks);


      lastOne->SetParameter(0,func->GetParameter(0));
      lastOne->SetParameter(1,func->GetParameter(1));
      lastOne->SetParameter(2,func->GetParameter(2));
      // lastOne->SetParLimits(2,0*func->GetParameter(2),0.5*func->GetParameter(2));

      lastOne->SetParameter(3,func->GetParameter(3));
      lastOne->SetParameter(4,func->GetParameter(4));
      lastOne->SetParameter(5,func->GetParameter(5));

      lastOne->SetParameter(6,func->GetParameter(6));
      lastOne->SetParameter(7,func->GetParameter(7));
      aux = 0;
      for(Int_t i = 1; i<n_peaks; i++){
        lastOne->SetParameter((i+8-1),func->GetParameter(i+8+aux));
        lastOne->SetParName((i+8-1),Form("A_{%d}",i+2));
        aux=aux+2;
      }

    
      lastOne->SetParName(0,"A_{baseline}");
      lastOne->SetParName(1,"#mu_{baseline}");
      lastOne->SetParName(2,"#sigma_{baseline}");
      lastOne->SetParName(3,"A_{1}");
      lastOne->SetParName(4,"#mu_{1}");
      lastOne->SetParName(5,"#sigma_{1}");
      lastOne->SetParName(6,"A_{2}");
      lastOne->SetParName(7,"#mu_{2}");
      // lastOne->SetParName(8,"#sigma_{2}");
    
      if(darknoise){
        lastOne->FixParameter(4,sphe_charge);
        lastOne->FixParameter(7,sphe_charge2);
      }
      if(fixZero){
        lastOne->FixParameter(0,0);
        lastOne->FixParameter(1,0);
        lastOne->FixParameter(2,1);
      }
      Int_t fit_status = -1;
      fit_status = hcharge->Fit("lastOne","R");
      cout << "Fit status before free std dev: " << fit_status << endl;

      if(make_free_stddevs == false){
        fu[0]->SetParameter(0,abs(lastOne->GetParameter(0)));
        fu[0]->SetParameter(1,lastOne->GetParameter(1));
        fu[0]->SetParameter(2,lastOne->GetParameter(2));
    
        fu[1]->SetParameter(0,abs(lastOne->GetParameter(3)));
        fu[1]->SetParameter(1,lastOne->GetParameter(4));
        fu[1]->SetParameter(2,lastOne->GetParameter(5));
    
        fu[2]->SetParameter(0,abs(lastOne->GetParameter(6)));
        fu[2]->SetParameter(1,lastOne->GetParameter(7));
        fu[2]->SetParameter(2,(TMath::Power((2),0.5)*lastOne->GetParameter(5)));
    
        for(Int_t i = 1; i<n_peaks; i++){
          fu[i+2]->SetParameter(0,abs(lastOne->GetParameter(i+7)));
          fu[i+2]->SetParameter(1,(lastOne->GetParameter(4) + (i+1)*(lastOne->GetParameter(7)-lastOne->GetParameter(4))));
          fu[i+2]->SetParameter(2,(TMath::Power((i+2),0.5)*lastOne->GetParameter(5)));
        }
      }
      if (debug_level == 3){
        lastOne->Draw("SAME");
        for(Int_t i = 0; i<(2+n_peaks); i++){
          fu[i]->SetLineColor(kGray+1);
          fu[i]->SetNpx(1000);
          fu[i]->Draw("SAME");
        }
        return;
      }

      MyFunctionFree MyFuncFree;
      MyFuncFree.n_peaks = n_peaks;

      TF1 *lastOneFree = new TF1("lastOneFree",MyFuncFree,xmin,xmax,7+n_peaks*2);

      if(make_free_stddevs == true){
        setParametersFree(lastOneFree,lastOne);
        // Int_t lastpar = 7+2*n_peaks-1;
        // Double_t lastGausUpperLim = 1.2*lastOneFree->GetParameter(lastpar);
        // Double_t lastGausLowerLim = 0.8*lastOneFree->GetParameter(lastpar);
        // lastOneFree->SetParLimits(lastpar,lastGausLowerLim,lastGausUpperLim);
        // cout << lastGausLowerLim << endl;
        // cout << lastGausUpperLim << endl;

        fit_status = hcharge->Fit("lastOneFree","R");
        cout << "Fit status with free std dev: " << fit_status << endl;


        fu[0]->SetParameter(0,lastOneFree->GetParameter(0));
        fu[0]->SetParameter(1,lastOneFree->GetParameter(1));
        fu[0]->SetParameter(2,lastOneFree->GetParameter(2));

        fu[1]->SetParameter(0,lastOneFree->GetParameter(3));
        fu[1]->SetParameter(1,lastOneFree->GetParameter(4));
        fu[1]->SetParameter(2,lastOneFree->GetParameter(5));

        fu[2]->SetParameter(0,lastOneFree->GetParameter(6));
        fu[2]->SetParameter(1,lastOneFree->GetParameter(7));
        fu[2]->SetParameter(2,lastOneFree->GetParameter(8));

        for(Int_t i = 1, j = 1; i<n_peaks; i++){
          fu[i+2]->SetParameter(0,abs(lastOneFree->GetParameter(j+8)));
          fu[i+2]->SetParameter(1,(lastOneFree->GetParameter(4) + (i+1)*(lastOneFree->GetParameter(7)-lastOneFree->GetParameter(4))));
          fu[i+2]->SetParameter(2,lastOneFree->GetParameter(j+1+8));
          j+=2;
        }

      }


      lastOne->SetNpx(1000);
      lastOneFree->SetNpx(1000);

      hcharge->Draw("hist");
      

      xmin = fu[0]->GetParameter(1)-5*fu[0]->GetParameter(2);

      hcharge->GetXaxis()->SetRangeUser(1.2*xmin,1.1*xmax);
      hcharge->StatOverflows(kTRUE);
      lastOne->SetRange(xmin,xmax);
      lastOneFree->SetRange(xmin,xmax);


      for(Int_t i = 0; i<(2+n_peaks); i++){
        fu[i]->SetLineColor(kGray+1);
        fu[i]->SetNpx(1000);
        fu[i]->Draw("SAME");
      }
      if(make_free_stddevs==false) lastOne->Draw("LP SAME");
      else {
        lastOneFree->Draw("LP SAME");
        lastOne = (TF1*)lastOneFree->Clone("lastOne");
      }
      string name = histogram + ".root";
      //     c1->Print(name.c_str());
      cout << "1th peak = " << lastOne->GetParameter(4) << endl;
      cout << "2th peak = " << lastOne->GetParameter(7) << endl;
      cout << "sphe charge = " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
      cout << " SNR = " << (lastOne->GetParameter(4))/sqrt(pow(lastOne->GetParameter(2),2)+pow(lastOne->GetParameter(5),2)) << endl;
      cout << " SNR2 = " << abs((lastOne->GetParameter(4))/lastOne->GetParameter(2)) << endl;
      out <<  lastOne->GetParameter(4) << " " << lastOne->GetParameter(7) << " " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;

      // ____________________________ Finish of sphe fit ____________________________ //


      //_________________ Drawing lines for the cross-talk probability _________________ //

      sphe_charge = lastOne->GetParameter(4);
      sphe_charge2 = lastOne->GetParameter(7);

      Double_t delta1 = (sphe_charge2 - sphe_charge)/deltaminus;
      Double_t delta2 = deltaplus*(sphe_charge2 - sphe_charge);
      //     Double_t delta2 = sphe_charge+(sphe_charge2 - sphe_charge)/2;

      Double_t ymax = hcharge->GetMaximum();

      TLine *l1 = new TLine(delta1,0,delta1,ymax);
      TLine *l2 = new TLine(delta2,0,delta2,ymax);
      l1->SetLineWidth(2);
      l2->SetLineWidth(2);
      l1->SetLineColor(kRed);
      l2->SetLineColor(kRed);
    
      l1->Draw("");
      l2->Draw("");
    
    
      if(darknoise){
        // ____________________________ Start dark noise CT analysis ____________________________ //
        
        cout << "\n\n\nMaking ratio between 1 sphe and 2 sphe: " << endl;
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram
        // could be taken, but the effort
        // is bigger.
        // Also, remember to multiply by
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);
        
        cout << "Total 1 = " << oneIntegral << " \t total 2 =  " << twoIntegral << endl;
        cout << "Ratio = " << twoIntegral/oneIntegral << endl;
        
        cout << "Total events normalized = " << lastOne->Integral(xmin,xmax) << endl;
        
        cout << "CT probability = " << twoIntegral/(oneIntegral+twoIntegral) << endl;
        
        
        
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        
        cout << "\n\n\n Another possibility would be CT = " << (twoIntegral+threeIntegral)/(oneIntegral+twoIntegral+threeIntegral) << endl;
        
        // ____________________________ Finish Poisson analysis ____________________________ //
        
      }
    
    
    
      if(is_poisson_test){
        
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram
        // could be taken, but the effort
        // is bigger.
        // Also, remember to multiply by
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        

        Double_t lambda = -TMath::Log(zeroIntegral/totalIntegral);
        
        cout << "Lambda = " << lambda << endl;
        
        Double_t ct = hcharge->GetMean()/(lambda*oneMean);
        
        cout << "Cross-talk = " << ct << endl;
        
        
        
        //         cout << "teste -> " << hcharge->GetMean() << endl;
        // ____________________________ Finish Poisson analysis ____________________________ //
      }
    }

    void setParametersFree(TF1 *lastOneFree, TF1 *lastOne){
      lastOneFree->SetParameter(0,lastOne->GetParameter(0));
      lastOneFree->SetParameter(1,lastOne->GetParameter(1));
      lastOneFree->SetParameter(2,lastOne->GetParameter(2));
      // lastOneFree->SetParLimits(2,0*lastOne->GetParameter(2),0.5*lastOne->GetParameter(2));

      lastOneFree->SetParameter(3,lastOne->GetParameter(3));
      lastOneFree->SetParameter(4,lastOne->GetParameter(4));
      lastOneFree->SetParameter(5,lastOne->GetParameter(5));

      lastOneFree->SetParameter(6,lastOne->GetParameter(6));
      lastOneFree->SetParameter(7,lastOne->GetParameter(7));
      lastOneFree->SetParameter(8,TMath::Power((2),0.5)*lastOne->GetParameter(5));

      for(Int_t i = 1, j= 1; i<n_peaks; i++){
        lastOneFree->SetParameter((j+9-1),lastOne->GetParameter(i+8-1));
        lastOneFree->SetParameter((j+1+9-1),(TMath::Power((i+2),0.5)*lastOne->GetParameter(5)));
        Double_t lastGausUpperLim = (1+stdVar)*lastOneFree->GetParameter(j+1+9-1);
        Double_t lastGausLowerLim = stdVar*lastOneFree->GetParameter(j+1+9-1);
        lastOneFree->SetParLimits(j+1+9-1,lastGausLowerLim,lastGausUpperLim);
        lastOneFree->SetParName((j+9-1),Form("A_{%d}",i+2));
        lastOneFree->SetParName((j+1+9-1),Form("#sigma_{%d}",i+2));
        // cout << "parameter " << j+9-1 << " = " << lastOneFree->GetParameter(j+9-1);
        // cout << " ... parameter " << j+1+9-1 << " = " << lastOneFree->GetParameter(j+1+9-1) << endl;
        j+=2;
      }


      lastOneFree->SetParName(0,"A_{baseline}");
      lastOneFree->SetParName(1,"#mu_{baseline}");
      lastOneFree->SetParName(2,"#sigma_{baseline}");
      lastOneFree->SetParName(3,"A_{1}");
      lastOneFree->SetParName(4,"#mu_{1}");
      lastOneFree->SetParName(5,"#sigma_{1}");
      lastOneFree->SetParName(6,"A_{2}");
      lastOneFree->SetParName(7,"#mu_{2}");
      lastOneFree->SetParName(8,"#sigma_{2}");

      if(darknoise){
        lastOneFree->FixParameter(4,sphe_charge);
        lastOneFree->FixParameter(7,sphe_charge2);
      }
      if(fixZero){
        lastOneFree->FixParameter(0,0);
        lastOneFree->FixParameter(1,0);
        lastOneFree->FixParameter(2,1);
      }
    }



    Calibration(string mname = "c"){
      myname = mname;
    }
    
};









// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //





class SPHE2{

  public:


    ///////////////////////////////////////////////////////////////////////////////
    //                      Variables to be changed by user                      //
    ///////////////////////////////////////////////////////////////////////////////

    Bool_t led_calibration = true; // if external trigger + led was used, set true
                                   // start and finish will be the time of integration

    Bool_t just_a_test = false; // well.. it is just a test, so `just_this` is the total waveforms analysed
    Int_t  just_this   = 200;
    Int_t  channel     = 1;
    string rootfile    = "analyzed.root";


    Double_t tolerance    = 5;  // n sigmas (smoothed) (not used for led)
                                // if `method` is set to `fix`, the threshold will be the absolute value of tolerance, no baseline is calculated
    Double_t baselineTime = 5000; // is the time to compute the baseline (not used for led)
                                  // If the `method` is set to `dynamic` the baseline is calculated over the range of baselineTime
                                  // and it is updated depending on the next point with a weigth of `influece`
                                  // If `method` is set to `static`, baseline is calculated once using baseLimit as cut
    Double_t baseLimit    = 3;  // higher then this wont contribute to the baseline abs(baseLimit) (not used for led)

    Double_t start  = 0;        // start the search for peaks or start the integration (led)
    Double_t finish = 10250;    // fisish the search or finish the integration (led)

    Double_t timeLow        = 60; // integration time before peak (not used for led)
    Double_t timeHigh       = 400; // integration time after peak (not used for led)
    Double_t lowerThreshold = -1; // threshold to detect noise (normal waveform) (not used for led)
    Double_t maxHits        = 1; // maximum hit before discarding   (not used for led)
    Double_t too_big        = 1000; // if there is a peak > "too_big" .. wait "waiting" ns for next peak
    Double_t waiting        = 1000;
    Double_t filter         = 14; // one dimentional denoise filter (0 equal no filder)
    Double_t interactions   = 45; // for moving avarage filter (not used on led)
    Int_t    ninter         = 2; // N times that moving average is applied

    Double_t dtime = 4.;        // time step in ns


    Double_t get_wave_form = false; // for getting spe waveforms
    Double_t mean_before   = 120; // time recorded before and after the peak found
    Double_t mean_after    = 1000;
    Double_t sphe_charge   = 1809.52; // charge of 1 and 2 p.e. (use fit_sphe.C)
    Double_t sphe_charge2  = 3425.95;
    Double_t sphe_std      = 500;

    // coeficients to surround gaussian of 1 spe.
    // Gain             = (sphe_charge2 - sphe_charge)
    // spe's get events where charge < Gain*deltaplus  and charge < Gain/deltaminus
    // If deltaminus is set to zero, sphe_std*deltaplus will be used instead
    // This value can be checked with fit_sphe.C
    Double_t deltaplus  = 1.3;
    Double_t deltaminus = 1.5;


    // Not so common to change
    Double_t social_distance = 2; // demands that there is a minimum distance of social_distance * timeHigh between 2 consecutive peaks found
    string   method          = "static"; // `dynamic` or `derivative` evaluation of the baseline
                               // `fix` will not evaluate baseline and use raw threshold
                               // See tolerance, baselineTime and baselineLimit above

    bool check_selection = true; // uses(or not) variable `selection` to discard wvfs
    Bool_t withfilter = true; // Integrate in the filtered waveform


    ///////////////////////////////////////////////////////////////////////////////
    //                         Finish variables for user                         //
    ///////////////////////////////////////////////////////////////////////////////


    // Here are more general ones

    ANALYZER *z      = nullptr; // Main analyser
    TFile    *fout   = nullptr; // File to save histogram
    TFile    *fwvf   = nullptr; // File to save waveforms
    TTree    *twvf   = nullptr; // tree of waveforms
    ADC_DATA *sample = nullptr; // For saving waveforms

    string myname; // name of the declared class
    Bool_t derivate = false; // control if user has set to derivate
    Int_t kch; // current channel used
    Int_t nshow = 100; // will show `nshow` debugging waveforms up
    Int_t nshow_start = 0; // starting from here
    Bool_t getstaticbase = false; // control if `static` method was set
    Double_t threshold = 0;

    Bool_t get_this_wvf = true;
    Bool_t get_this_charge = true;

    TGraph *g_smooth   = nullptr; // for the nshow sample plots
    TGraph *g_normal   = nullptr;
    TGraph *g_selected = nullptr;
    TGraph *g_discarded = nullptr;
    Double_t *timeg = new Double_t[memorydepth];


    // ____________________ Variables to calculate and reset ____________________ //

    TH1D *hbase        = nullptr; // not used now
    TH1D *hbase_smooth = nullptr; // compute smoothed baseline
    TH1D *hcharge      = nullptr; // final spe histogram

    Double_t mean = 0; //mean and stddev for hbase
    Double_t stddev = 0;

    vector<Double_t> peakPosition; //store position of the peaks
    vector<Double_t> peakMax; // store peak maximum
    vector<Int_t> peaksFound; // peaks found (threshold apply)
    vector<Int_t> peaksRise; // Rising derivative
    vector<Int_t> peaksCross; // Falling derivative

    vector<Double_t> smooted_wvf; // not needed to reset this
    vector<Double_t> denoise_wvf;

    vector<Double_t> selected_peaks; // these are for plotting
    vector<Double_t> selected_time;
    vector<Double_t> selected_charge;
    vector<Double_t> selected_charge_time;


    vector<Double_t> discardedTime; //store position of thed
    vector<Double_t> discardedPeak; //store position of thed
    vector<Int_t> discarded_idx; //store idx (0, 1 or 2) of the peaks
    // ____________________ ________________________________ ____________________ //

    vector<string> discard_reason = {"Social distance", "Negative hits", "Too big"};

    SPHE2(string m_name) : myname{m_name}{
      hbase        = new TH1D(Form("hbase_sphe_%s",myname.c_str()),"histogram for baseline",5*800,-400,400);
      hbase_smooth = new TH1D(Form("hbase_smooth_sphe_%s",myname.c_str()),"histogram for baseline smoothed",5*800,-400,400);
      hcharge      = new TH1D(Form("hcharge_sphe_%s",myname.c_str()),Form("hcharge_sphe_%s",myname.c_str()),50000,0,0);

      smooted_wvf.resize(memorydepth);
      denoise_wvf.resize(memorydepth);


    }


    void searchForPeaks(){
    }

    /**
     * This function searches and integrate peaks, used for sphe
     * If `led_calibration` is set to true, integration is made between `start` and `finish`
     * Load everything
     * Create TTree/Histograms necessary
     * Reset area, clean vectors, histograms, etc. for next waveform
     * Process data: moving average, denoise filters, etc.
     * Search for peaks
     * Integrate peaks and store the waveform
     * Plot graphs for debugging
     * Save histogram
     */
    void giveMeSphe(){
      gROOT->SetBatch(kTRUE);

      if(!z){ // load analyzer in case it is nullptr
        z = new ANALYZER(myname.c_str());
        z->dtime = dtime;
        z->setAnalyzer(rootfile);
      }
      if(!z->t1){ // otherwise, check if the ttree was loaded
        z->dtime = dtime;
        z->setAnalyzer(rootfile);
      }

      z->kch = channel;
      kch = channel;
      fout = new TFile(Form("sphe_histograms_darkCount_Ch%i.root",kch),"RECREATE");

      // ____________________ Setting up what needs to be set ____________________ //

      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",kch),"RECREATE");
      }
      else{ // in the case we are not taking waveform, I change this values for the same setup of integration
        mean_before = timeLow;
        mean_after = timeHigh;
      }
      twvf = new TTree("t1","mean waveforms");
      // if the user set mean_before = 24 and mea50n_after = 100
      // we need to get  6 + 25 points + 1 = 32  (+1 because of the peak start)
      // if I found a point in around 100 ns, that is i = 25, I perform a while that goes
      // from i-mean_before = 19 up to i+mean_after = 50 (included) = 32 pts
      sample = new ADC_DATA();
      Int_t npts_wvf = (mean_before/dtime + mean_after/dtime) + 1;
      sample->setBranchName(npts_wvf);
      twvf->Branch(Form("Ch%i",channel),&sample,sample->tobranch.c_str());

      if(led_calibration==false){
        method = "led";
      }
      else{

        for(int i = 0; i < memorydepth; i++){
          timeg[i] = i*dtime;
        }
        peaksFound.reserve(memorydepth);
        // this might be too much, time will tell. Keep an eye in the RAM memory
        selected_peaks.reserve(memorydepth);
        selected_time.reserve(memorydepth);
        selected_charge.reserve(memorydepth);
        selected_charge_time.reserve(memorydepth);
      }

      if(method == "static") getstaticbase = true;
      else if(method == "fix"){
      }
      else if(method == "derivative"){
        derivate = true;
        peaksRise.reserve(memorydepth);
        peaksCross.reserve(memorydepth);
      }

      Double_t delta = sphe_charge2 - sphe_charge;
      if(deltaminus == 0){
        deltaminus = sphe_charge - deltaplus*sphe_std;
        deltaplus = sphe_charge + deltaplus*sphe_std;
      }
      else{
        deltaminus = delta/deltaminus;
        deltaplus = delta*deltaplus;
      }



      // ____________________ ___________________________ ____________________ //



      Int_t nentries = z->nentries;
      if(just_a_test){nentries = just_this;}
      for(Int_t i = 0; i<nentries; i++){
        theGreatReset();
        z->getWaveform(i,kch); // get waveform by memory
        if(check_selection && z->ch[kch].selection != 0) continue;

        processData();

        if(method == "static"){
          searchForPeaks();
        }
        else if(derivate){
          // Search for crossing points in the derivated waveform
          // positive crossing stored in peaksRise
          // negative crossing stored in peaksCross
          // from start to finish, discounting what was for getting integral with an extra point (not waveform, this filter is done later)
          z->zeroCrossSearch(&smooted_wvf[0], peaksRise, peaksCross, start+abs(timeLow)+dtime, finish-abs(mean_after)-dtime);
          derivateApplyThreshold();
        }
        cleanPeaks();
        integrateSignals();
        if(snap()){
          drawMySamples();
          cout << "Event " << z->ch[kch].event << " total of peaks: " << peaksFound.size() << ", Valid = " << selected_charge.size() << endl;
          for(unsigned int j = 0; j < selected_charge.size(); j++){
            cout << "\t\t Charge = " << selected_charge[j] << " at " << selected_charge_time[j] << endl;
          }
        } // draw sample graphs
      }


      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",kch));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
      }


    }
    void theGreatReset(){
        hbase->Reset();
        if(method != "led"){
          if(derivate){
            peaksRise.clear();
            peaksCross.clear();
          }
          else{
            hbase_smooth->Reset();
          }
          peakPosition.clear();
          peakMax.clear();
          selected_peaks.clear();
          selected_time.clear();
          selected_charge.clear();
          selected_charge_time.clear();
          discardedTime.clear();
          discardedPeak.clear();
          discarded_idx.clear();
          peaksFound.clear();
        }
      
        denoise_wvf.clear();
        smooted_wvf.clear();
    }

    void processData(){
      z->applyDenoise(filter, z->ch[kch].wvf, &denoise_wvf[0]); // denoise is stored at denoise_wvf. If no filter, they are equal

      /**
       * From https://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
       *It makes no difference whether the smooth operation is applied before or after the differentiation. What is important, however,
       *is the nature of the smooth, its smooth ratio (ratio of the smooth width to the width of the original peak), and the number of
       *times the signal is smoothed. The optimum values of smooth ratio for derivative signals is approximately 0.5 to 1.0. For a first
       *derivative, two applications of a simple rectangular smooth (or one application of a triangular smooth) is adequate. For a second
       *derivative, three applications of a simple rectangular smooth or two applications of a triangular smooth is adequate. The general rule is:
       *for the nth derivative, use at least n+1 applications of a rectangular smooth. (The Matlab signal processing program iSignal automatically
       *provides the desired type of smooth for each derivative order).
       *Usediscarded_idx to check the best option
       **/
      if(derivate){
        z->differenciate(1e3,z->ch[kch].wvf,&smooted_wvf[0]); // multiply by 1e3 so we can see something :)
      }
      for(Int_t i = 0; i < interactions; i++){
        if (i == interactions-1 && method == "static") getstaticbase = true;
        smooted_wvf = movingAverage(&smooted_wvf[0], interactions, getstaticbase);
      }

    }

    bool goodSocialDistance(Int_t id2, Int_t id1){
      if(abs(id2-id1)<= social_distance*timeHigh){
        return false;
      }
      else{
        return true;
      }
    }
    bool checkTooBig(Bool_t &wait_now, Double_t &refWait, Int_t pos){
      Double_t peakmax;
      if(derivate){
        peakmax = denoise_wvf[pos];
      }
      else{
        peakmax = peakMax[pos];
      }
      if(peakmax > too_big){
        wait_now = true;
        refWait = pos*dtime;
        return true;
      }
      else{
        wait_now = false;
        return false;
      }

    }

    void derivateApplyThreshold(){
      Int_t ntotal = (int)peaksCross.size();
      if(peaksRise.size() < ntotal){ // there should be only pairs
        cout << "@@@@@@@@@@@@@@@@@@@@ THIS SHOULD NOT HAPPEN @@@@@@@@@@@@@@@@@@@@" << endl;
        cout << "\n\n\n\n\n\n\n\n\n\n" << endl;
        cout << "@@@@@@@@@@@@@@@@@@@@ THIS SHOULD NOT HAPPEN @@@@@@@@@@@@@@@@@@@@" << endl;
        ntotal = peaksRise.size();
      }

      for(Int_t i = 0; i < ntotal; i++){
        Double_t crossPositive = peaksRise[i];
        Double_t candidatePosition = peaksCross[i]; // peaks previously selected
        if(crossPositive*dtime > social_distance*timeHigh){
          crossPositive = social_distance*timeHigh/dtime;
        }
        z->getMaximum(crossPositive*dtime, candidatePosition*dtime);
        if(z->temp_max < tolerance) {
          continue;
        }
        peaksFound.push_back(peaksCross[i]) ;
      }

    }
    void cleanPeaks(){
      Int_t ntotal = (int)peaksFound.size();
      Bool_t wait_now = false;
      Double_t refWait = 0;

      for(Int_t i = 0; i < ntotal; i++){
        Double_t candidatePosition = peaksFound[i]; // peaks previously selected
        if(wait_now){ // I am waiting for a big pulse to finish
          if(candidatePosition*dtime >= (refWait+waiting)){
            wait_now = false;
          }
          else{
            if(snap()){
              discardedTime.push_back(candidatePosition*dtime);
              discardedPeak.push_back(smooted_wvf[candidatePosition]);
              discarded_idx.push_back(2);
            }
            continue;
          }
        }

        //request social distance
        if(!goodSocialDistance(peaksFound[i+1],candidatePosition)){
          if(snap()){ // discarding current peak
            discardedTime.push_back(candidatePosition*dtime);
            discardedPeak.push_back(smooted_wvf[candidatePosition]);
            discarded_idx.push_back(0);
          }

          // discarding  next peak
          // But first! Check if it is not too big =3
          if(checkTooBig(wait_now, refWait, candidatePosition)){
            if(snap()){
              discardedTime.push_back(candidatePosition*dtime);
              discardedPeak.push_back(smooted_wvf[candidatePosition]);
              discarded_idx.push_back(2);
            }
            i++;
            continue;
          }
          if(snap()){
            discardedTime.push_back(peaksFound[i+1]*dtime);
            discardedPeak.push_back(smooted_wvf[peaksFound[i+1]]);
            discarded_idx.push_back(0);
          }
          i++;
          continue;
        } // and before
        else if(i != 0 && !goodSocialDistance(peaksFound[i-1], candidatePosition)){
          if(snap()){
            discardedTime.push_back(candidatePosition*dtime);
            discardedPeak.push_back(smooted_wvf[candidatePosition]);
            discarded_idx.push_back(0);
          }
          continue;
        }
        // social distance all good
        // Now, checking for big pulses
        if(checkTooBig(wait_now, refWait, candidatePosition)){
          if(snap()){
            discardedTime.push_back(candidatePosition*dtime);
            discardedPeak.push_back(smooted_wvf[candidatePosition]);
            discarded_idx.push_back(2);
          }
          continue;
        }

        peakPosition.push_back(peaksFound[i]);
      }
    }

    vector<Double_t> delay_line(vector<Double_t> v, Double_t delay_time){
      if(delay_time==0) return v;
      vector<Double_t> res(v.size());
      for(int i=0; i<(int)v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
      }
      return res;
    }

    template<class T>
    vector<Double_t> movingAverage(T* v, Int_t myinte, Bool_t eval_baseline){
      Int_t midpoint = 0;
      Int_t interactions = 0;
      Double_t width = 0;
      Double_t sum = 0;
      Int_t n = memorydepth;

      vector<Double_t> res(n,0);
      if(myinte==0) { // nothing to do
        for (Int_t i = 0; i < n; i++) {
          res[i] = v[i];
        }
        return res;
      }

      if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
      }
      else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
      }

      for(Int_t i = 0; i < n; i++){
        if(i<midpoint || i>(n-midpoint) || interactions == 0){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          res[i] = v[i];
        }
        else if(i*dtime > start && i*dtime < finish){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];
          }
         res[i] = (sum/width);

         if(eval_baseline){ // in case we are computing the baseline
           if(i*dtime<=baselineTime && abs(res[i])<baseLimit){
             hbase_smooth->Fill(res[i]);
           }
         }
        }
        else{
          res[i] = 0;
        }
        sum=0;
      }
      if(eval_baseline){ // in case we are computing the baseline
        mean = hbase_smooth->GetMean();
        stddev = hbase_smooth->GetStdDev();
      }
      return res;
    }

    void integrateSignals(){
      unsigned int npeaks = peakPosition.size();
      for(unsigned int i = 0; i < npeaks; i++){
        Int_t peakIdx = 0;
        if(!led_calibration) peakIdx = peakPosition[i];
        getIntegral(peakIdx); // If get_mean_wvf is set to false, there is no problem!
        Double_t charge = z->temp_charge;
        Double_t peak = z->temp_max;
        if(get_this_charge){
          hcharge->Fill(charge);
          if(snap())
          {
            selected_charge.push_back(charge);
            selected_charge_time.push_back(peakIdx*dtime);
          }
        }
        else{
          continue; // If I didn't get the charge, I will not take the waveform
        }
        if(get_this_wvf && get_this_wvf){
          if(charge>=deltaminus && charge <= deltaplus){
            sample->selection = 1;
          }
          sample->peak = peak;
          sample->charge = charge;
          twvf->Fill();
        }
      }
    }

    void getIntegral(Int_t peakIdx){
      get_this_charge = true;
      get_this_wvf = true;
      Double_t from = peakIdx - mean_before/dtime;
      Double_t to = peakIdx + mean_after/dtime;
      if(from < 0 || to >= memorydepth){ // check if boundaries are being respected
        from = peakIdx - timeLow/dtime;
        to = peakIdx + timeHigh/dtime;
        if(get_wave_form) get_this_wvf = false; // in the case not, I dont take the waveform
      }
      Int_t iwvf = 0;

      Double_t res = 0;
      Double_t max = -1e12;
      Int_t negativeHits = 0;
      Double_t integralfrom = peakIdx - timeLow/dtime;
      Double_t integralto = peakIdx - timeHigh/dtime;

      if(led_calibration){
        from = 0;
        to = memorydepth-1;
        integralfrom = start/dtime;
        integralto = finish/dtime;
        get_this_wvf = true;
      }
      for(Int_t i = from; i <= to; i++, iwvf++){
      Double_t val = denoise_wvf[i];
        if(get_wave_form && get_this_wvf){
          sample->wvf[iwvf] = z->ch[kch].wvf[i];
        }
        if(i >= integralfrom && i <= integralto){
          if (val < lowerThreshold){
            negativeHits++;
          }
          if(negativeHits >= maxHits){
            get_this_wvf = false;
            get_this_charge = false;
            if(snap()){
              discardedTime.push_back(peakIdx*dtime);
              discardedPeak.push_back(smooted_wvf[peakIdx]);
              discarded_idx.push_back(1);
            }
            break;
          }
          if(snap())
          {
            selected_peaks.push_back(val);
            selected_time.push_back(i*dtime);
          }
          if(withfilter == false){ val = z->ch[kch].wvf[i]; }
          res += val;
          if(val>=max){
            max = val;
          }
        }
      }
      // store values here because of lazyness
      z->temp_charge = res*dtime;
      z->temp_max = max;

    }

    // ____________________________________________________________________________________________________ //
    Bool_t snap(){
      if(z->currentEvent >= nshow_start && z->currentEvent < nshow_start + nshow && !led_calibration){
        return true;
      }
      else{
        return false;
      }
      return false; // just to be sure oO
    }
    void drawMySamples(){
      string sampleName = "ev_" + to_string(z->ch[kch].event)+"_"+to_string(kch);

      g_smooth = new TGraph(memorydepth, timeg, &smooted_wvf[0]);
      g_normal = new TGraph(memorydepth, timeg, &denoise_wvf[0]);
      // TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str(),1920,0,700,500);
      // this is not working when saving
      TCanvas *c1 = new TCanvas();

      c1->cd(1);
      g_smooth->SetLineColor(kRed);
      g_smooth->SetLineWidth(3);

      g_normal->SetLineColor(kBlue);
      g_normal->SetTitle(" ");

      g_normal->Draw("AL");
      g_smooth->Draw("L SAME");

      if (derivate) {
        threshold = tolerance;
        mean = 0;
      }

      TLine *lmean = new TLine(start,mean,memorydepth*dtime,mean);
      TLine *ldev = new TLine(finish,mean+threshold,memorydepth*dtime,mean+threshold);

      lmean->SetLineColor(kGreen);
      ldev->SetLineColor(kGreen);

      lmean->SetLineWidth(2);
      ldev->SetLineWidth(2);

      lmean->Draw();
      ldev->Draw();
      Int_t n = selected_peaks.size();
      if(n!=0){
        g_selected = new TGraph(n, &selected_time[0], &selected_peaks[0]);
        g_selected->SetMarkerColor(kBlack);
        g_selected->Draw("P* SAME");
      }

      Int_t ndis = discardedPeak.size();
      if(ndis != 0){
        g_discarded = new TGraph(ndis, &discardedTime[0], &discardedPeak[0]);
        g_discarded->SetMarkerColor(kRed);
        g_discarded->Draw("P* SAME");
      }



      // fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");


    }



};























class SPHE{
       
    
    
    
  public:
    
    TFile *fout=nullptr;

    TTree *tout = nullptr;

    Bool_t creation = true; //to verify creation of ttree branches


    TFile *fwvf = nullptr;
    TTree *twvf = nullptr;


    Double_t value = 0;
    Double_t desv = 0;
    Int_t channel = 1;
    ADC_DATA ch;
    ADC_DATA sample;


    // ____________________ Variables to calculate ____________________ //

    TH1D *hbase = new TH1D("hbase","histogram for baseline",5*800,-400,400);
    TH1D *hbase_smooth = new TH1D("hbase_smooth","histogram for baseline smoothed",5*800,-400,400);
    // TH1D *hcharge = new TH1D("hcharge","",100000,-50000,50000);
    TH1D *hcharge = new TH1D("hcharge","",50000,0,0);
    TH1D *hzero = new TH1D("hzero","",120000,-200000,2*1300000);
    TH1D *hnobase = new TH1D("hnobase","",120000,-200000,2*1300000);

    TH1D *hstat = new TH1D("hstat","",120000,-200000,2*1300000);


    Double_t dtime = 2.;

    Double_t mean = 0;
    Double_t stddev = 0;
    Double_t tolerance; // n sigmas
    Double_t threshold;
    Bool_t fix_threshold = false;
    Double_t baseLimit;

    Double_t filter = 0;
    Bool_t withfilter = true;
  
    Double_t timeLimit; // time after LED signal
    Double_t timeLow ; // integration time before peak
    Double_t timeHigh ; // integration time after peak

    Double_t start = 0;
    Double_t finish = memorydepth*dtime;

    Double_t noiseLow;
    Double_t noiseHigh;

    Double_t forBaseline = 2000;
    Double_t baseCharge = 0;
    
    Double_t social_distance = 2;
    Bool_t check_selection = true;

    Int_t interactions; // for moving avarage
    Int_t midpoint = 0; // midpoint that takes the avarage
    Int_t width = 0;
    Double_t sum=0;

    vector<Double_t> peakPosition;
    vector<Double_t> peakMax;

    Double_t charge = 0;
    Double_t treeCharge = 0;
    Double_t treePeak = 0;
    Double_t strikes = 0;
    Double_t ptsLow = 0;
    Double_t ptsHigh = 0;
    Bool_t discard = false;

    TH1D *hfilter = new TH1D("hfilter","",5*800,-400,400);
    Double_t devFilter;


    Int_t cut = 0;

    Int_t gap = 2000; // at the and, we wont look at the final 4000 ns
    Double_t lowerThreshold = -9999;
    Double_t maxHits = 15;
    Double_t nSigmas_integration = 0.5;
    Int_t atleast = 120; //at least this time in the first integration and second integration
    Int_t lowPtsLimit = 20;
    Bool_t just_a_test = false; //if true, will only look 5000 events;
    Bool_t led_calibration = false; // for led integration only

    Int_t just_this = 5000;

    // ____________________ matching filter ____________________

    vector<Double_t> mySample;
    vector<Double_t> matched;
    vector<Double_t> matchingAreas;

    Double_t shift = 0;
    Double_t matched_value = 0;
    Double_t higherValue = 0;
    Double_t tempPosition;


    Double_t area_off = 0;
    Double_t matchTrigger = 1.;

    Double_t from = 0;
    Double_t to = 0;


    Bool_t matching = false;
    Bool_t badBaseline = false; // this need to be used in the first Gamma run...

    // ____________________ Variables to show events ____________________ //

    const static Int_t nshow = 100;
    Double_t my_events[nshow];
    Double_t eventNow = 0;
    Int_t aux_events = 0;


    TGraph *g_smooth = nullptr;
    TGraph *g_normal = nullptr;
    TGraph *g_points = nullptr;
    
    vector<Double_t> temp_peak;
    vector<Double_t> peak_smooth;
    vector<Double_t> timeg;

    vector<Double_t> selected_peaks;
    vector<Double_t> selected_time;
    Int_t teste = 0;


    Bool_t fillg = true;

    Int_t mark = 0;

    Double_t xi = 0; // time interval for sampling (ns)
    Double_t xf = 1000;
    Double_t center = 300;
    Int_t nsample = (xf-xi)/2+1;
    Double_t counter = 0;

    TH2D *hsample = new TH2D("sample","peak vs time", nsample, xi, xf, 3000,-1500,1500);

    vector<Double_t> ymean;
    vector<Double_t> ynormal;
    vector<Double_t> xmean;

    // variables to find mean waveform
    Bool_t get_wave_form = false;
    vector<Double_t> mean_waveforms;
    Double_t wvfcharge;
    Bool_t valid;
    Int_t naverages;
    Double_t mean_before = 100;
    Double_t mean_after = 400;
    Int_t low_cut = 800; // to select sphe waveforms, inside this region, peak should not be bigger then val_cut
    //it was 920 before
    Int_t high_cut = 2000;
    Double_t val_cut = 6;
    TMultiGraph *gm = new TMultiGraph();
    TGraph *gwaveforms = nullptr;
    Double_t shifter = 12;

    Double_t sphe_charge;
    Double_t sphe_charge2;
    Double_t delta;
    Double_t deltaplus = 1.4;
    Double_t deltaminus = 1.2;
    Double_t sphe_std;

    Double_t sphe_charge_ch0;

    Double_t sphe_charge2_ch0;

    Double_t sphe_std_ch0;
    Double_t sphe_std_ch1;

    string charge_status = "";

    Bool_t darkNoise = false;

    Double_t baselineTime=2000;
    Double_t too_big = 1000;
    Double_t waiting = 0;
    Double_t pre_filter = 0;
    Int_t peakPerEvent = 0;
    // ____________________________________________________________________________________________________ //

    void giveMeSphe_darkCount(string name){
      temp_peak.resize(memorydepth);
      gROOT->SetBatch(kTRUE);
      fout = new TFile(Form("sphe_histograms_darkCount_Ch%i.root",channel),"RECREATE");
      tout = new TTree("t1","baseline info");
      darkNoise = true;
      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
      }
      twvf = new TTree("t1","mean waveforms");

      if(led_calibration==false){
        makeHistogram(name);
        fout->WriteObject(tout,"t1","TObject::kOverwrite");
        // fout->Close();
      }
      else{
        makeSimpleHistogram(name);
        // get_wave_form = false;
      }
      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",channel));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
      }
      channel=channel + 1;
    }


    // ____________________________________________________________________________________________________ //
    void giveMeSphe(string name){
      gROOT->SetBatch(kTRUE);
  
      fout = new TFile("sphe_histograms.root","RECREATE");
      tout = new TTree("t1","baseline info");
      darkNoise = false;
      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
      }
      twvf = new TTree("t1","mean waveforms");
  
      makeHistogram(name);

      fout->WriteObject(tout,"t1","TObject::kOverwrite");
      fout->Close();
      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",channel));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
      }
      channel=channel + 1;
    }

    // ____________________________________________________________________________________________________ //
    void integrateSignal(){
      charge = 0;
      Int_t npeaks = (int)peakPosition.size();
      value = 0;
      desv = 0;

      Int_t discardedPeaks = 0;

      Int_t aux_sample = 0;
    

      sphe_charge = sphe_charge_ch0; // wave0
      sphe_charge2 = sphe_charge2_ch0; // wave0
      delta = sphe_charge2 - sphe_charge;
      sphe_std = sphe_std_ch0;

      Double_t statcharge[npeaks];
      Double_t statpeak[npeaks];
      Double_t trackLow[npeaks];
      Double_t trackHigh[npeaks];
      Double_t trackStrikes[npeaks];
      Double_t statpos[npeaks];
      Bool_t discard_this[npeaks];
    
      vector<vector<Double_t>> temp_waveforms(npeaks);
      vector<Double_t> waveforms(memorydepth);
      //     vector<TGraph*> gwaveforms(npeaks);
      vector<Bool_t> notAGoodWaveform(npeaks);


      Int_t auxstat = 0;
      Bool_t notGood = false;
      Double_t lastTime = 0 ;
      Double_t peakStd = 0;
    
      Double_t thismanypointsLow = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t thismanypointsHigh = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t thismanypointsBase = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t totalmany = 0; //this helps speacially for matching...
      Bool_t next_is_bad = false;
      Double_t waitingInterval = -1;
      for(Int_t i = 0; i<npeaks; i++){
        statcharge[i] = 0;
        statpeak[i] = 0;
        trackStrikes[i] = 0;
        discard_this[i] = false;
        notGood = false;
        thismanypointsLow = 0;
        thismanypointsHigh = 0;
        thismanypointsBase = 0;
        totalmany = 0;
      
      
        if((peakPosition[i]+social_distance*timeHigh>=peakPosition[i+1] && i+1!=npeaks) || next_is_bad || peakPosition[i]<waitingInterval){

          selected_peaks.push_back(temp_peak.at(peakPosition[i]/dtime));
          selected_time.push_back(peakPosition[i]);
          discard_this[i] = true;
          for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
            if(temp_peak.at(j)>=statpeak[i]){
              statpeak[i] = temp_peak.at(j);
            }
          }
          if(statpeak[i]>too_big){
         
            waitingInterval = peakPosition[i]+waiting;
          }
          if(peakPosition[i]<waitingInterval && waitingInterval!=-1) next_is_bad = true;
          else next_is_bad = false;
          continue;
        }
        else{
          waitingInterval = -1;
          next_is_bad = false;
        }
      
      
      
        notAGoodWaveform[i] = false;

            
                
        Int_t strikes = 0;
        // integration for the back part
        for(Int_t j = (peakPosition.at(i))/dtime; j>= (peakPosition[i]-timeLow)/dtime; j--){
          charge += temp_peak.at(j);
          if(withfilter) statcharge[i]+= temp_peak.at(j);
          else statcharge[i]+= ch.wvf[j];
          if(temp_peak.at(j)>=statpeak[i]){
            statpeak[i] = temp_peak.at(j);
          }
          selected_peaks.push_back(temp_peak.at(j));
          selected_time.push_back(timeg.at(j));
          thismanypointsLow++;
          totalmany++;
          if(temp_peak[j]<lowerThreshold){
            strikes++;
          }
        }
      
      
        // integration for the front part
        for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
          charge += temp_peak.at(j);
        
          if(withfilter) statcharge[i]+= temp_peak.at(j);
          else statcharge[i]+= ch.wvf[j];
          if(temp_peak.at(j)>=statpeak[i]){
            statpeak[i] = temp_peak.at(j);
          }
          selected_peaks.push_back(temp_peak.at(j));
          selected_time.push_back(timeg.at(j));
          thismanypointsHigh++;
          totalmany++;
          if(temp_peak[j]<lowerThreshold){
            strikes++;
          }
        }
      
        if(statpeak[i]>too_big){
          waitingInterval = peakPosition[i]+waiting;
          discard_this[i] = true;
        }
        if(strikes>=maxHits){
          //         cout << "strikes " << strikes << ".... " << eventNow << endl;
          discard_this[i] = true;
        }
      
        // for get_wave_form
        if(peakPosition.at(i) - mean_before>=0 && peakPosition.at(i)+mean_after<memorydepth*dtime){
          for(Int_t j = peakPosition.at(i)/dtime - mean_before/dtime; j <= peakPosition.at(i)/dtime+mean_after/dtime; j++){
            // temp_waveforms[i].push_back(temp_peak.at(j));
            temp_waveforms[i].push_back(ch.wvf[j]);
            //           if(temp_waveforms[i].size()>200/dtime && temp_waveforms[i].size()<400/dtime){
            //             if(temp_peak.at(j)<-20) notAGoodWaveform[i]=true;
            //           }
          }
        }
        else{
          notAGoodWaveform[i] = true;
        }
      
        trackHigh[i] = thismanypointsHigh;
        trackLow[i] = thismanypointsLow;
      
      
      
      
        if(snap()){
        
          cout << "charge = " << statcharge[i]*dtime << " at " << peakPosition.at(i) << " with " << thismanypointsLow << " + " << thismanypointsHigh << " or this " << thismanypointsBase << " discard ? " << discard_this[i] << endl;
        }
      
      
      
        
        
        
      }
      for(Int_t i = 0; i<npeaks; i++){
        if(discard_this[i]==false){

          peakPerEvent++;
          statcharge[i] = dtime*statcharge[i];
          charge = dtime*charge;
          hcharge->Fill(statcharge[i]);
          hnobase->Fill(statcharge[i]);
          hstat->Fill(statcharge[i]);
          treeCharge = statcharge[i];
          treePeak = statpeak[i];
          ptsHigh = trackHigh[i];
          ptsLow = trackLow[i];
          charge_status += " charge = " + to_string(statcharge[i]);
          tout->Fill();
          //             if (charge<-800){
          //               cout << "\n\n------>event = " << eventNow << " " << charge << " " << peakPosition[i] << endl;
          //             }
          //             if(statcharge[i]>10000 && statcharge[i]<11000 && snap()){
          //                 cout << ".............................. charge = " << statcharge[i] << " at " << peakPosition.at(i) << " event = " << eventNow << endl;
          //             }
            
          if(get_wave_form){
            Double_t newbase = 0;
            if(notAGoodWaveform[i]==false){
              if(statcharge[i]>=delta/deltaminus && statcharge[i]<=delta*deltaplus){
                naverages++;
                valid = true;
                high_cut = mean_before+mean_after;
                for(Int_t j = low_cut/dtime; j<high_cut/dtime; j++){
                  if(temp_waveforms[i][j]>val_cut){
                    naverages--;
                    valid = false;
                    break;
                  }
                }
              }
              else{
                valid = false;
              }
              for(Int_t j = 0; j<(int)waveforms.size(); j++){
                if(j<(int)temp_waveforms[i].size()){
                    
                  // waveforms[j] = temp_waveforms[i][j];
                  waveforms[j] = temp_waveforms[i][j];
                  sample.wvf[j]=waveforms[j];
                  if(valid) mean_waveforms[j]+=waveforms[j];
                }
                else{
                  waveforms[j] = 0;
                  sample.wvf[j]=0;
                  mean_waveforms[j]+=waveforms[j];
                }
                  
              }
              wvfcharge = statcharge[i];
              sample.event = eventNow;
              sample.charge = statcharge[i];
              sample.selection = valid;
              twvf->Fill();
                
              //gwaveforms[i] = new TGraph(waveforms.size(),&timeg[0],&waveforms[0]);
              //gm->Add(gwaveforms[i],"LP");
            }
          }
            
            
        }
        else{
          //             charge_status += " Discarded - > " + to_string(statcharge[i]*dtime) + " ";
          //             charge_status += " strikes - > " + to_string(trackStrikes[i]) + " ";
        }
      }
    
    

    

 
    }

    // ____________________________________________________________________________________________________ //
    void searchForPeaks(){
      // For each peak found, we i am looking for the maximum above tolerance
      Int_t n = (int)peak_smooth.size();
      Int_t npeaks = 0;
   
      higherValue = 0;
      tempPosition = 0;
      Bool_t rebase = false;
      Bool_t normalbase = true; // make it true to calculate
      Double_t gapstart = baselineTime;
      Double_t gapend = 6000;
      Double_t actual_finish = finish;
      //     for(Int_t i = gapstart/dtime; i<=gapend/dtime; i++){
      //       if(temp_peak[i]>200){
      //         actual_finish = gapstart;
      //       }
      //     }

      threshold = tolerance*stddev;
      if(fix_threshold) threshold = tolerance;
      if(check_selection && ch.selection!=0) return;
                              
      for(Int_t i = 0; i<n; i++){
      
        if(i>timeLow/dtime && i<(n-timeHigh/dtime) && (i>=start/dtime && i<=actual_finish/dtime)){ //searching out of the beginning
          // if(i>baselineTime/dtime && temp_peak[i]>300){ // not sure why this was here
          // if(i>baselineTime/dtime){
          //   break;
          // }
          if(peak_smooth[i]>(mean+threshold) && peak_smooth[i-1]<(mean+threshold) && (i<=gapstart || i>=gapend)){
            npeaks++;
          
            peakMax.push_back(peak_smooth.at(i));
            peakPosition.push_back(timeg.at(i));

            if(snap()){
              cout << "npeaks = " << npeaks << " " << peakPosition.at(npeaks-1) << " " << peakMax.at(npeaks-1) << " " << eventNow << endl;
            
            }
            //           if(i<=baselineTime/dtime){
            //               rebase = true;
            //           }
            //           if(rebase && normalbase){
            //             rebase =false;
            //             normalbase=false;
            //             Int_t binmax = hbase->GetMaximumBin();
            //             Double_t newB1 = hbase->GetXaxis()->GetBinCenter(binmax);
            //             Double_t newB2 = hbase->GetMean();
            //             Double_t newB = (newB1<newB2)? newB1 : newB2;
            //
            //             for(Int_t k = 0; k<temp_peak.size();k++){
            //               temp_peak[k] = temp_peak[k]-newB;
            // //               mean = newB;
            //               if(eventNow==my_events[nshow-1]){
            // //                 cout << eventNow << " " << k << " " << temp_peak[k] << " " << newB << " " << hbase->GetMean() << endl;
            //               }
            //             }
            //           }
          
          }
        
        
  
        }
        else if(i>finish/dtime){
          break;
        }
      }
    
    }

    vector<Double_t> delay_line(vector<Double_t> v, Double_t delay_time){
      if(delay_time==0) return v;
      vector<Double_t> res(v.size());
      for(int i=0; i<(int)v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
      }
      return res;
    }

    // ____________________________________________________________________________________________________ //
    void lookWaveform(){

      peakPerEvent = 0;
    
      vector<Double_t> ma_to_shift = movingAverage(&temp_peak[0],pre_filter);
      vector<Double_t> shifted=delay_line(ma_to_shift, shifter);//cusp(MIN[i], h);
      smoothWithMovingAvarage(shifted); // first calculate avarage

 
      searchForPeaks(); //search for peaks with the moving avarage
      if(snap()){ // if the events match, create the graphs
      
        g_smooth = new TGraph(timeg.size(),&timeg[0],&peak_smooth[0]);
        g_normal = new TGraph(timeg.size(),&timeg[0],&temp_peak[0]);
        //         cout << mean << " " << stddev << endl;
      
      }
      integrateSignal();
    
      if(snap()){
        drawMySamples();
        cout << "Event " << eventNow << " total of peaks: " << peakPerEvent << endl;
      } // draw sample graphs
    }

    // ____________________________________________________________________________________________________ //
    void makeHistogram(string filename){


      if(matching == true){
        getMySample();
        
        if(shift==0){
          shift = (to-from)/2; //make sure to not use "2." we want an integer here
        }
        if(area_off==0){area_off = to - from;}
        
        cout << area_off << " " << shift << endl;
      }

    
      if(creation){
        tout->Branch("value",&value,"value/D");
        tout->Branch("desv",&desv,"desv/D");
        tout->Branch("channel",&channel,"channel/D");
        tout->Branch("treeCharge",&treeCharge,"treeCharge/D");
        tout->Branch("treePeak",&treePeak,"treePeak/D");
        tout->Branch("ptsLow",&ptsLow,"ptsLow/D");
        tout->Branch("ptsHigh",&ptsHigh,"ptsHigh/D");
        tout->Branch("strikes",&strikes,"strikes/D");
        //         creation = false;

        twvf->Branch(Form("Ch%i",channel),&sample,sample.tobranch.c_str());
      }
    
      cout << "reading: " << filename << endl;
      string rootfile = filename + ".root";
    
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");

      TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
      bch->SetAddress(&ch);
      Int_t nentries = t1->GetEntries();
    
      TH1D *hdark = new TH1D("hdark","",1000,-10000,30000);

      //_____________________________ Start creating random data to check _____________________________ //
      TRandom *rmd = new TRandom();

      my_events[0] = 0;
      for(Int_t i = 1; i<nshow; i++){
        //         my_events[i] = static_cast<Int_t>(rmd->Uniform(my_events[i-1]+1,my_events[i-1]+50));
        my_events[i] = i;
      }
      sort(my_events,my_events+nshow);
    


      eventNow = 0;
      aux_events = 0;
      //_____________________________ End creating random data to check _____________________________ //
    
    
    
      // ____________________ Start resetting globals ____________________ //
      mean = 0;
      stddev = 0;
      midpoint = 0; // midpoint that takes the avarage
      width = 0;
      sum=0;
      charge = 0;
      discard = false;
      eventNow = 0;
      aux_events = 0;
      teste = 0;
      fillg = true;
      hsample->Reset();
      counter = 0;
      // ____________________ Finish resetting globals ____________________ //
    
    

    
      Double_t aux = 0;
      Int_t j = 0;
    
      temp_peak.clear();
      temp_peak.resize(memorydepth);
      timeg.clear();
      peak_smooth.clear();
      peakPosition.clear();
      peakMax.clear();
      selected_peaks.clear();
      selected_time.clear();
    
    
      matched.clear();
    
      ymean.resize(nsample);
      ynormal.resize(nsample);
      xmean.resize(nsample);
      for(Int_t i = 0; i<(int)ymean.size(); i++){
        ymean.at(i) = 0;
        xmean.at(i) = 0;
      }

    
      hbase->Reset();
      hstat->Reset();
      hbase_smooth->Reset();
      hcharge->Reset();
      hzero->Reset();
      hnobase->Reset();
    
      mean_waveforms.clear();
      mean_waveforms.resize(memorydepth,0);
      naverages = 0;
    
      for(Int_t i = 0; i<memorydepth; i++){
        timeg.push_back(dtime*i);
      }
    
      if(just_a_test){nentries = just_this;}
      //     aux = 0;
   
      for(Int_t i = 0; i<nentries; i++){
        bch->GetEvent(i);
        DENOISE dn;
        if(filter>0)dn.TV1D_denoise<Double_t>(&ch.wvf[0],&temp_peak[0],memorydepth,filter);
        else{
          for(Int_t i = 0; i<memorydepth; i++){
            temp_peak[i] = ch.wvf[i];
          }
        }
        // for(Int_t i = 0; i<memorydepth; i++){
        //   if((i<=5000/dtime) && abs(temp_peak[i])<6){
        //     hbase->Fill(temp_peak[i]);
        //   }
        // }
      
      
        aux = ch.event;
        if(static_cast<Int_t>(eventNow)%200==0){
          cout << eventNow << "\r" << flush;
        }
        eventNow =  ch.event;
        // here is were all the calculation is done !!!
        lookWaveform();
      
        if(snap() && aux_events+1<nshow){
          aux_events = aux_events+1;
        }
      
        hbase->Reset();
        hbase_smooth->Reset();
      
        charge_status = "";
        temp_peak.clear();
        temp_peak.resize(memorydepth);
        peak_smooth.clear();
        peakPosition.clear();
        peakMax.clear();
        selected_peaks.clear();
        selected_time.clear();
        baseCharge = 0;
        matched.clear();
      
      }
    
    
      cout << "____________________________ " << counter/nsample << " \n \n \n" << endl;
    
      if(get_wave_form){
        for(Int_t i = 0; i<memorydepth; i++){
          //             cout << i << " " << mean_waveforms[i] << " ";
             
          mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
          //             cout << mean_waveforms[i] << endl;
            
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;
        
      }
    
      TCanvas *cwvf = new TCanvas();
      cwvf->cd();
      //     gm->Draw("A");
    
      TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);
    
      if(get_wave_form){
        //         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
      }

      fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

      f1->Close();

    }

    // ____________________________________________________________________________________________________ //
    void smoothWithMovingAvarage(vector<Double_t> shifted){
    
      Int_t n = shifted.size();
      if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
      }
      else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
      }
    
    
      for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-midpoint) || interactions == 0){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          peak_smooth.push_back(shifted.at(i));
        }
        else if(i>cut/dtime){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+shifted.at(j);
            //                 cout << sum << endl;
          }
          peak_smooth.push_back(sum/width);
            
          if(timeg.at(i)<=baselineTime && abs(peak_smooth.at(i))<baseLimit){
            hbase_smooth->Fill(peak_smooth.at(i));
          }

        }
        else{
          peak_smooth.push_back(0);
        }

        
        
        sum=0;
      }
      mean = hbase_smooth->GetMean();
      stddev = hbase_smooth->GetStdDev();
    
    
    
    }
    // ____________________________________________________________________________________________________ //
    void drawMySamples(){

      string sampleName = "ev_" + to_string(static_cast<Int_t>(my_events[aux_events]))+"_"+to_string(channel);
  
      // TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str(),1920,0,700,500);
      // this is not working when saving
      TCanvas *c1 = new TCanvas();

      c1->cd(1);
      g_smooth->SetLineColor(kRed);
      g_smooth->SetLineWidth(3);
    
      g_normal->SetLineColor(kBlue);
      g_normal->SetTitle(charge_status.c_str());
    
      g_normal->Draw("AL");
      g_smooth->Draw("L SAME");

      TLine *lmean = new TLine(timeLimit,mean,memorydepth*dtime,mean);
      TLine *ldev = new TLine(timeLimit,mean+threshold,memorydepth*dtime,mean+threshold);
    
      lmean->SetLineColor(kGreen);
      ldev->SetLineColor(kGreen);
    
      lmean->SetLineWidth(2);
      ldev->SetLineWidth(2);
    
      lmean->Draw();
      ldev->Draw();
      Int_t n = selected_peaks.size();
      if(n!=0){
        g_points = new TGraph(n,&selected_time[0],&selected_peaks[0]);
        g_points->SetMarkerColor(kBlack);
        g_points->Draw("P* SAME");
      }
    
    
    
    
      // Some fancy drawing now;
      if(matchingAreas.size()!=0){
        Int_t nAreas = matchingAreas.size(); //always two points per area of course
        TLine *larea;
        Double_t max = g_normal->GetYaxis()->GetXmax();
        for(Int_t i = 0; i<nAreas; i++){
          larea = new TLine(matchingAreas[i],0,matchingAreas[i],max);
          cout << "\t\t\t" << matchingAreas[i] << " ";
          if((i+1)%2==0)cout << "\n" << endl;
          larea->SetLineColor(kRed);
          larea->SetLineWidth(2);
          larea->Draw();
        }
      }
    
      fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");

    
    }



    // ____________________________________________________________________________________________________ //
    Bool_t snap(){
      if(eventNow == my_events[aux_events]){
        return true;
      }
      else{
        return false;
      }
      return false; // just to be sure oO
    }

    void getMySample(){
      TFile *fsample = new TFile("mysample.root","READ");
      if (fsample->IsZombie()) {
        matching = false;
        std::cout << "\n\n\n\n\n Error opening file \n\n\n\n\n" << std::endl;
        return;
      }
      TTree *tsample = (TTree*)fsample->Get("t1");
    
      Double_t values;
      tsample->SetBranchAddress("values",&values);
      Int_t nsample = tsample->GetEntries();
      for(Int_t i = 0; i < nsample; i++){
        tsample->GetEntry(i);
        mySample.push_back(values);
      }
      fsample->Close();
      return;
    
    }

    Bool_t checkAreas(Double_t totalmany){ //return true if the points did not touched the area...
      Int_t nAreas = matchingAreas.size(); //always two points per area of course
      Int_t n = selected_time.size();
      for(Int_t i = n; i>n-totalmany; i--){
        for(Int_t j = 0; j<nAreas; j++){
          if(selected_time[i-1]>=matchingAreas[j] && selected_time[i-1]<=matchingAreas[j+1]){
            return true;
          }
          j++;
        }
      }
      return false;
    }

    void makeSimpleHistogram(string filename){


      sphe_charge = sphe_charge_ch0; // wave0
      sphe_charge2 = sphe_charge2_ch0; // wave0
      delta = sphe_charge2 - sphe_charge;
      sphe_std = sphe_std_ch0;
  
      if(creation){
        twvf->Branch(Form("Ch%i",channel),&sample,sample.tobranch.c_str());
      }
      cout << "reading: " << filename << endl;
      string rootfile = filename + ".root";
  
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");

      TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
      bch->SetAddress(&ch);
      Int_t nentries = t1->GetEntries();
      Double_t charge = 0;
      Bool_t noise = false;
      Int_t noise_hits = 0;
      Double_t max = -1e12;
      mean_waveforms.clear();
      mean_waveforms.resize(memorydepth,0);
      naverages = 0;
    
    
      ofstream ftmp;
      ftmp.open("valid_events.log",ios::out);
      if(just_a_test){nentries = just_this;}
      for(Int_t i = 0; i<nentries; i++){
        bch->GetEvent(i);
        DENOISE dn;
        dn.TV1D_denoise<Double_t>(&ch.wvf[0],&temp_peak[0],memorydepth,filter);

        noise = false;
        max = -1e12;
        for(Int_t j = start/dtime; j<finish/dtime; j++){
          if(withfilter) charge += temp_peak.at(j);
          else charge += ch.wvf[j];
          if(temp_peak[j]>=max){
            max = temp_peak[j];
          }
          if(temp_peak[j]>too_big) noise=true;
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }
        for(Int_t j = 0; j<start/dtime-800/4; j++){
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }
        for(Int_t j = finish/dtime+800/4; j<memorydepth; j++){
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }

    
    
        if(noise==false){
          if(check_selection && ch.selection!=0){
          }
          else{            
            valid = false;
            for(Int_t j = 0; j<memorydepth; j++){
              sample.wvf[j] = ch.wvf[j];
            }
            if(charge*dtime>=delta/deltaminus  && charge*dtime<=delta*deltaplus){
              valid = true;
              for(Int_t j = 0; j<memorydepth; j++){
                mean_waveforms[j]+=ch.wvf[j];
              }
              naverages++;
            }

            wvfcharge = charge*dtime;
            sample.charge = wvfcharge;
            sample.selection = valid;
            twvf->Fill();
            hcharge->Fill(charge*dtime);
            // cout << charge*dtime << endl;
            ftmp << i << "\n";
          }
        }
        charge=0;
      }
      ftmp.close();

      if(get_wave_form){
        for(Int_t i = 0; i<memorydepth; i++){
          timeg.push_back(dtime*i);
          mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;
        
      }

  
      TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);
      if(get_wave_form){
        //         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
      }
      fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

      f1->Close();
    
    }

    template<class T>
    vector<Double_t> movingAverage(T* v, Int_t myinte){
    
      Int_t n = memorydepth;
      vector<Double_t> res(n,0);
      if(myinte==0) {
        for (Int_t i = 0; i < n; i++) {
          res[i] = v[i];
        }
        return res;
      }
      if(myinte%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = myinte/2+1;    //midpoint will be 5 here
        width = myinte+1;
      }
      else{
        midpoint = (myinte-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = myinte;
      }
    

      for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          res[i] = v[i];
        }
        else if(i>cut/dtime){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];
            //                 cout << sum << endl;
          }
          res[i] = (sum/width);
            
        
        }
        else{
          v[i] = (0);
        }
        
        sum=0;
      }
      return res;
    
    }
  


};




class MeanSignal{
  
  public:
  
    Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits
  
    vector<Int_t> channels = {1,2};
    Double_t minval = 0;
    Double_t maxval = 1e12;
    Double_t avoid_saturation = 1e6;
    Double_t fprompt = 0.8; // well... fprompt
    string mustbe = "bigger"; // must be 'bigger' or 'lower' then fprompt set.
  
    Bool_t checkFprompt(Double_t fprompt_data, Double_t fprompt_ref, string mustbe){
      if(mustbe == "bigger"){
        if(fprompt_data>=fprompt_ref) return true;
        else return false;
      }
      else{
        //     cout << "Test ... " << fprompt_data << " " << fprompt_ref << endl;
        if(fprompt_data<=fprompt_ref){
          return true;
        }
        else return false;
      }
    }
    void mean_signal(string filename){
  
      string rootfile = filename + ".root";
  
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");
      vector<ADC_DATA> ch(channels.size());
      vector<TBranch*> bch(channels.size());
      for(Int_t k = 0; k<(int)channels.size();k++){
        bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
        bch[k]->SetAddress(&ch[k]);
      }
      Int_t nentries = t1->GetEntries();
      vector<Double_t> norm(channels.size(),0);
  
      vector<vector<Double_t>> avg(channels.size(),std::vector<Double_t>(memorydepth,0));
      vector<vector<Double_t>> avgn(channels.size(),std::vector<Double_t>(memorydepth,0));
      vector<Double_t> time(memorydepth);
      for(Int_t j = 0; j<memorydepth; j++){
        time[j]+=j*dtime;
      }
  
      TFile *fout = new TFile("averaged_waveforms.root","RECREATE");
      vector<TH2D*> hpersistence(channels.size());
      for(Int_t k = 0; k<(int)channels.size();k++) hpersistence[k] = new TH2D(Form("hpersistence_%i",channels[k]),Form("hpersistence_%d",channels[k]),5000,0,20000,500,-500,avoid_saturation);

      for(Int_t i = 0; i<nentries; i++){
        for(Int_t k = 0; k<(int)channels.size();k++){
          bch[k]->GetEvent(i);
      
          if(ch[k].charge>minval && ch[k].charge<maxval && ch[k].peak<avoid_saturation){
            if(checkFprompt(ch[k].fprompt,fprompt,mustbe)){
              norm[k]+=1;
              for(Int_t j = 0; j<memorydepth; j++){
                hpersistence[k]->Fill(j*dtime,ch[k].wvf[j]);
                avg[k][j]+=ch[k].wvf[j];
              }
            }
          }
        }
      }
  
      vector<TGraph*> gavg(channels.size());
      vector<TGraph*> gavgn(channels.size());
      vector<Double_t> maxvalue(channels.size());
      for(Int_t k = 0; k<(int)channels.size(); k++){
        maxvalue[k] = *std::max_element(begin(avg[k]),end(avg[k]))/norm[k];
        for(Int_t j = 0; j<memorydepth; j++){
          avg[k][j]=avg[k][j]/norm[k];
          avgn[k][j]=avg[k][j]/maxvalue[k];
        }
        cout << "Ch" << channels[k] << " total waveforms = " << norm[k] << endl;
        gavg[k] = new TGraph(memorydepth,&time[0],&avg[k][0]);
        gavgn[k] = new TGraph(memorydepth,&time[0],&avgn[k][0]);
        fout->WriteObject(gavg[k],Form("average_ch%i",channels[k]),"TObject::kOverwrite");
        fout->WriteObject(gavgn[k],Form("average_normalized_ch%i",channels[k]),"TObject::kOverwrite");
        fout->WriteObject(hpersistence[k],Form("hpersistence_%i",channels[k]),"TObject::kOverwrite");
      }
  
  
  
    }
  
};
















class Resolution{
  
  public:

    vector<Int_t> channels = {1,2};
    TFile *fout = nullptr;
    TTree *tout = nullptr;
    Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits

    Double_t intmin = 5600;
    Double_t intmax = 20000;
    Double_t offset = 400; // just so it is easier
    Double_t intstep = 1000;
    Int_t nints = (intmax-(intmin+offset))/intstep+4;
    vector<Double_t> integrations;
    vector<vector<Double_t>> resCharge;
    vector<Double_t> lowlim = {450,530,630,700,730,820,850,850,850,850,850,850,850,850,850,850,850,850};
    vector<Double_t> average = {500,550,600,800,850,850,900,950,950,950,950,950,950,950,950,950,950,950};
    Double_t ampl = 7e4;
    void trackResolution(string filename){
      integrations.resize(nints);
      resCharge.resize(channels.size(),std::vector<Double_t>(nints,0));
  
      string rootfile = filename + ".root";
  
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");
      vector<ADC_DATA> ch(channels.size());
      vector<TBranch *> bch(channels.size());
      for(Int_t k = 0; k<(int)channels.size(); k++){
        bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
        bch[k]->SetAddress(&ch[k]);
      }
      Int_t nentries = t1->GetEntries();
  
      fout = new TFile("resolution.root","RECREATE");
      tout = new TTree("resolution","tout");
  

      integrations[0] = offset+100;
      integrations[1] = integrations[0]+100;
      integrations[2] = integrations[1]+200;
      integrations[3] = integrations[2]+200;
      integrations[4] = offset+intstep;
      for(Int_t i = 5; i<nints; i++) integrations[i] = integrations[i-1]+intstep;
      vector<TH1D*> hspecs(nints);
      vector<TF1*> fa(nints);


      for(Int_t i = 0; i<nints; i++){
        hspecs[i] = new TH1D(Form("spec_%.0f",integrations[i]),Form("t = %.0f ns",integrations[i]),400,-100,1900);
        fa[i] = new TF1(Form("fa_run%.0f",integrations[i]),"([0]/(2*[2]))*exp((x-[1])/[2]+[3]*[3]/(2*[2]*[2]))*TMath::Erfc(((x-[1])/[3]+[3]/[2])/TMath::Power(2.,0.5))",lowlim[i],1300);
      }
  
  
  
      ifstream fsphe;
      fsphe.open("sphe.txt",ios::in);
      Double_t temp;
      vector<Double_t> sphes(channels.size());
      for(Int_t k = 0; k<(int)channels.size(); k++){
        fsphe >> temp >> temp >> temp;
        sphes[k] = temp;
      }
      fsphe.close();
  
  
      Double_t charge1,charge2;
      Double_t photoelec;
      for(Int_t i = 0; i<nentries; i++){
        //   for(Int_t i = 0; i<1; i++){
        cout << "reading event " << i << "\r" << flush;
        photoelec = 0;
        for(Int_t k = 0; k<(int)channels.size(); k++){
          bch[k]->GetEvent(i);
          integrate(ch[k].wvf,resCharge[k]);
        }
        if(ch[0].fprompt>0.5 && ch[1].fprompt>0.5 && ch[0].peak<2800 && ch[1].peak<2800){
          for(Int_t j = 0; j<nints; j++){
            for(Int_t k = 0; k<(int)channels.size(); k++){
              photoelec += resCharge[k][j]/sphes[k];
            }
            hspecs[j]->Fill(photoelec);
            photoelec=0;
          }
        }


    
      }
      f1->Close();
  
      vector<Double_t> sigmu(nints);
      vector<Double_t> ersigmu(nints);
      vector<Double_t> mu(nints);
      vector<Double_t> ermu(nints);
      THStack *hs = new THStack();
      for(Int_t j = 0; j<nints; j++){
        fa[j]->SetParameters(ampl,average[j],150,30);
        fa[j]->SetParName(0,"A");
        fa[j]->SetParName(1,"#mu");
        fa[j]->SetParName(2,"#tau");
        fa[j]->SetParName(3,"#sigma");
        hspecs[j]->Fit(Form("fa_run%.0f",integrations[j]),"R0Q");
        hspecs[j]->Fit(Form("fa_run%.0f",integrations[j]),"R0Q");
        fa[j]->SetRange(0,1300);
        sigmu[j] = fa[j]->GetParameter(3)/fa[j]->GetParameter(1);
        ersigmu[j] = sigmu[j]*sqrt(pow(fa[j]->GetParError(3)/fa[j]->GetParameter(3),2)+pow(fa[j]->GetParError(1)/fa[j]->GetParameter(1),2));
        mu[j] = fa[j]->GetParameter(1);
        ermu[j] = fa[j]->GetParError(1);
        fout->WriteObject(hspecs[j],Form("spec_%.0f",integrations[j]),"TObject::kOverwrite");
        fout->WriteObject(fa[j],Form("spectrum_run%.0f",integrations[j]),"TObject::kOverwrite");
        hspecs[j]->SetLineWidth(2);
        hs->Add(hspecs[j],"hist");
    
      }
  
      //   gStyle->SetPalette(kVisibleSpectrum);
      TGraphErrors *g = new TGraphErrors(nints,&integrations[0],&sigmu[0],0,&ersigmu[0]);
      g->GetXaxis()->SetTitle("Integration time (ns)");
      g->GetYaxis()->SetTitle("#sigma/#mu");
      g->SetMarkerStyle(21);
      g->SetMarkerSize(0.7);
      TCanvas *c1 = new TCanvas();
      g->Draw("AP");
      fout->WriteObject(g,"resolution","TObject::kOverwrite");
      TGraphErrors *g2 = new TGraphErrors(nints,&integrations[0],&mu[0],0,&ermu[0]);
      g2->GetXaxis()->SetTitle("Integration time (ns)");
      g2->GetYaxis()->SetTitle("#mu");
      g2->SetMarkerStyle(21);
      g2->SetMarkerSize(0.7);
      TCanvas *c2 = new TCanvas();
      g2->Draw("AP");
      fout->WriteObject(g2,"peaks","TObject::kOverwrite");
      TCanvas *c3 = new TCanvas();
      hs->Draw("nostack plc");
      c3->BuildLegend();
    
    }

    void integrate(Double_t wvf[], vector<Double_t> &charges){
      Double_t refCharge = 0;
      Int_t index_ints = 0;
      for(Int_t i = 0; i<memorydepth; i++){
        if(i>=intmin/dtime && i<(intmin+integrations[index_ints])/dtime){
          refCharge+=wvf[i]*dtime;
          if(i==memorydepth-1){
            charges[index_ints] = refCharge;
            index_ints++;
          }
      
        }
        else if(i==(intmin+integrations[index_ints])/dtime){
          //       cout << "filling i = " << i*dtime << " "  << integrations[index_ints] << " with " << refCharge << endl;
          charges[index_ints] = refCharge;
          refCharge+=wvf[i]*dtime;
          index_ints++;
        }
    
      }
  
    }
  
  
  
};









































class TimeDistribuction{
  public:
  
    vector<Int_t> channels={1,2};
    Double_t dtime = 4;
    Int_t interactions = 15;
    Bool_t just_a_test = false;
    Int_t just_this = 200;
    Int_t nshow = 20;
    Double_t lowerThreshold = 10;
    Double_t gap = 48; // in ns
    vector<Double_t> gtime;
  
  
    vector<double> delay_line(Double_t v[], Int_t delay_time){
      vector<double> res(memorydepth);
      for(int i=0; i<memorydepth; i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
      }
      return res;
    }

    void drawGraphs(TFile *f, TGraph * g1, TGraph * g2, Int_t count,Int_t k){
      TCanvas *c1 = new TCanvas();
      g1->SetLineColor(kBlue);
      g2->SetLineColor(kRed);
      g1->Draw("ALP");
      g2->Draw("LP SAME");
      f->WriteObject(c1,Form("g%i_ch%i",count,channels[k]),"TObject::kOverwrite");
    }
  
    
    // ____________________________________________________________________________________________________ //
    void time_distribuction(string filename){
      gROOT->SetBatch(kTRUE);
      Int_t nchannels = channels.size();
      vector<vector<Double_t>> peaks(channels.size(),std::vector<Double_t>(memorydepth,0));
      vector<vector<Double_t>> shifted(channels.size(),std::vector<Double_t>(memorydepth,0));
      vector<TGraph *> gsample(nshow*nchannels);
      vector<TGraph *> gshifted(nshow*nchannels);
      vector<TH1D*> htd(nchannels);
      vector<Int_t> posaux(nchannels);
      vector<vector<Double_t>> position(nchannels);
  
      string rootfile = filename + ".root";
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");
  
      TFile *fout = new TFile("time_distribuction.root","RECREATE");
  
      vector<ADC_DATA>  ch(nchannels);
      vector<TBranch *> bch(nchannels);
      for(Int_t k = 0; k<nchannels;k++){
        bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
        bch[k]->SetAddress(&ch[k]);
        htd[k] = new TH1D(Form("h_ch%i",channels[k]),Form("h_ch%i",channels[k]),memorydepth,0,dtime*memorydepth);
      }
      for(Int_t j = 0; j<memorydepth; j++){
        gtime.push_back(dtime*j);
      }
  
  
      Int_t nentries = t1->GetEntries();
      if(just_a_test) nentries = just_this;
      for(Int_t i = 0; i<nentries; i++){
        //     cout << "scanning event: " << i  << endl;
        cout << "scanning event: " << i << "\r" << flush;
        for(Int_t k = 0; k<nchannels;k++){
          bch[k]->GetEvent(i);
        }
    
    
        for(Int_t k = 0; k<nchannels;k++){
          shifted[k]=delay_line(ch[k].wvf, 40);
          smoothWithMovingAvarage(shifted[k]);
        }
        for(Int_t k = 0; k<nchannels;k++){
          for(Int_t j = 1; j<memorydepth-1; j++){
            if(shifted[k][j]>=lowerThreshold ){
              if(shifted[k][j]>shifted[k][j-1]){
                position[k].push_back(searchMax(ch[k].wvf,shifted[k],j)*dtime);
                posaux[k]++;
              }
            }
          }
        }
        if(i<nshow){
          for(Int_t k = 0; k<nchannels;k++){
            gsample[nchannels*i+k] = new TGraph(memorydepth,&gtime[0],&ch[k].wvf[0]);
            gshifted[nchannels*i+k] = new TGraph(memorydepth,&gtime[0],&shifted[k][0]);
          }
        }
      }
  
      cout << "\n";
      Int_t aux = 0;
  
      for(Int_t k = 0; k<(int)channels.size(); k++){
        for(Int_t i = 0; i<(int)position[k].size();i++){
          htd[k]->Fill(position[k][i]);
        }
        fout->WriteObject(htd[k],Form("h_ch%i",channels[k]),"TObject::kOverwrite");
      }

      for(Int_t k = 0; k<(int)channels.size(); k++){
        aux=0;
        for(Int_t i = 0; i<nshow; i++){
          //       cout << "scanning event: " << channels.size()*i+k << "\r" << flush;
          drawGraphs(fout,gsample[channels.size()*i+k],gshifted[channels.size()*i+k],aux,k);
          aux++;
        }
      }
      gROOT->SetBatch(kFALSE);
  
  
    }

    Int_t searchMax(Double_t wvf[],vector<Double_t> ref,Int_t &j){
      Double_t res = -1e12;
      Int_t maxpos = 0;
      for(Int_t i = j; ;i++){
        if(i<memorydepth){
          if(wvf[i]>res){
            res = wvf[i];
            maxpos = i;
          }
        }
        if(ref[i]<=lowerThreshold){
          j = i;
          break;
        }
      }
      return maxpos;
    }
  
  
    // ____________________________________________________________________________________________________ //
    void smoothWithMovingAvarage(vector<Double_t> &shifted){
      vector<Double_t> peak_smooth;
      Int_t n = (int)shifted.size();
      Double_t sum = 0;
      Int_t midpoint;
      Double_t width;
      if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
      }
      else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
      }
  
      for(Int_t i = 0; i < n; i++){
    
        if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          peak_smooth.push_back(shifted.at(i));
        }
        else if(i>0){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+shifted.at(j);
            //                 cout << sum << endl;
          }
          peak_smooth.push_back(sum/width);
        }
        else{
          peak_smooth.push_back(0);
        }
        sum=0;
      }
  
      for(Int_t i = 0; i<memorydepth; i++){
        shifted[i] = peak_smooth[i];
      }
  
  
    }
  
  
  
  
  
  
  
  
};




















































































// backup of previous if

// if(statcharge[i]>=delta/2. && statcharge[i]<=delta*1.5 && notAGoodWaveform[i]==false){
//   naverages++;
//   for(Int_t j = 0; j<(int)waveforms.size(); j++){
//     if(j<temp_waveforms[i].size()){
//       if(j<=50){
//         newbase += temp_waveforms[i][j];
//       }
//       else if(j==51){
//         
//         newbase/=51;
//         for(Int_t k=0; k<=51; k++){
//           waveforms[k] = temp_waveforms[i][k]-newbase;
//           wvf[k]=waveforms[k];
//           mean_waveforms[k]+=waveforms[k];
//         }
//       }
//       else{
//         waveforms[j] = temp_waveforms[i][j]-newbase;
//         wvf[j]=waveforms[j];
//         mean_waveforms[j]+=waveforms[j];
//       }
//     }
//     else{
//       waveforms[j] = 0;
//       wvf[j]=0;
//       mean_waveforms[j]+=waveforms[j];
//     }
//     
//   }
//   twvf->Fill();
//   
//   gwaveforms[i] = new TGraph(waveforms.size(),&timeg[0],&waveforms[0]);
//   gm->Add(gwaveforms[i],"LP");
// }
              
