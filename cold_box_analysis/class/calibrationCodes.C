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

    TH1D *htemp = nullptr;
    TH1D *hcharge = nullptr;
  
  
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
        h = (TH1D*)htemp->Clone("");
      }


      h->Rebin(rebin);
      Double_t scale = 1/(h->Integral());
      // h->Scale(scale);
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
        h->GetYaxis()->SetTitle("# of entries");
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
      TFile *f1 = nullptr;
      if (rootFile != "") {
        TFile *f1 = new TFile(rootFile.c_str(),"READ");
        hcharge = (TH1D*)f1->Get(histogram.c_str());
      }
      else{
        hcharge = (TH1D*)htemp->Clone("");
      }
      hcharge->Rebin(rebin);
    
      ofstream out;
      out.open("sphe.txt", ios::app);

      // ____________________________ Start of sphe fit ____________________________ //
      // hcharge->GetYaxis()->SetTitle("Normalized count");
      hcharge->GetYaxis()->SetTitle("# of events");
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
      //hcharge->Scale(scale);
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


        fu[0]->SetParameter(0,abs(lastOneFree->GetParameter(0)));
        fu[0]->SetParameter(1,lastOneFree->GetParameter(1));
        fu[0]->SetParameter(2,lastOneFree->GetParameter(2));

        fu[1]->SetParameter(0,abs(lastOneFree->GetParameter(3)));
        fu[1]->SetParameter(1,lastOneFree->GetParameter(4));
        fu[1]->SetParameter(2,lastOneFree->GetParameter(5));

        fu[2]->SetParameter(0,abs(lastOneFree->GetParameter(6)));
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
      cout << "Std dev  = " << lastOne->GetParameter(5) << endl;
      cout << "sphe charge = " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
      cout << " SNR = " << (lastOne->GetParameter(4))/sqrt(pow(lastOne->GetParameter(2),2)+pow(lastOne->GetParameter(5),2)) << endl;
      cout << " SNR2 = " << abs((lastOne->GetParameter(4))/lastOne->GetParameter(2)) << endl;
      out <<  lastOne->GetParameter(4) << " " << lastOne->GetParameter(7) << " " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;

      // ____________________________ Finish of sphe fit ____________________________ //


      //_________________ Drawing lines for the cross-talk probability _________________ //

      sphe_charge = lastOne->GetParameter(4);
      sphe_charge2 = lastOne->GetParameter(7);
      Double_t sphe_std = lastOne->GetParameter(5);

      Double_t delta1;
      Double_t delta2;
      if(deltaminus!=0)
      {
        delta1 = (sphe_charge2 - sphe_charge)/deltaminus;
        delta2 = deltaplus*(sphe_charge2 - sphe_charge);
      }
      else{
        delta1 = (sphe_charge-deltaplus*sphe_std);
        delta2 = (sphe_charge+deltaplus*sphe_std);
      }
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

    Bool_t led_calibration = false; // if external trigger + led was used, set true
                                   // start and finish will be the time of integration

    Bool_t just_a_test = false; // well.. it is just a test, so `just_this` is the total waveforms analysed
    Int_t  just_this   = 200;
    Int_t  channel     = 1;
    string filename    = "analyzed";

    vector<Int_t> nshow_range = {0,100}; // will save some debugging waveforms inside the range.
                                // Example, nshow_range = {0,2} will save waveforms 0, 1 and 2;

    Double_t tolerance    = 5;  // n sigmas (smoothed) (not used for led)
                                // if `method` is set to `fix`, the threshold will be the absolute value of tolerance, no baseline is calculated
    Double_t baselineTime = 5000; // is the time to compute the baseline (not used for led)
                                  // If the `method` is set to `dynamic` the baseline is calculated over the range of baselineTime
                                  // and it is updated depending on the next point with a weigth of `influece`
                                  // If `method` is set to `static`, baseline is calculated once using baseLimit as cut
    Double_t baseLimit    = 3;  // higher then this wont contribute to the baseline abs(baseLimit) (not used for led)

    Double_t start  = 0;        // start the search for peaks or start the integration (led)
    Double_t finish = 10250;    // fisish the search or finish the integration (led)

    Double_t timeLow         = 60; // integration time before peak (not used for led)
    Double_t timeHigh        = 400; // integration time after peak (not used for led)
    Double_t lowerThreshold  = -1; // threshold to detect noise (normal waveform) (not used for led)
    Double_t maxHits         = 1; // maximum hit before discarding   (not used for led)
    Double_t too_big         = 1000; // if there is a peak > "too_big" .. wait "waiting" ns for next peak
    Double_t waiting         = 1000;
    Double_t filter          = 14; // one dimentional denoise filter (0 equal no filder)
    Double_t interactions    = 45; // for moving avarage filter (not used on led)
    Int_t    ninter          = 2; // N times that moving average is applied
    Double_t diff_multiplier = 1e3; //derivative can be very small. Use this to make it easier to see
    Bool_t derivate_raw      = true; // If you apply a denoise in the data and want to take derivative on that, set to false

    Double_t dtime = 4.;        // time step in ns


    Double_t get_wave_form = false; // for getting spe waveforms
    Double_t mean_before   = 120; // time recorded before and after the peak found
    Double_t mean_after    = 1000;
    Double_t sphe_charge   = 1809.52; // charge of 1 and 2 p.e. (use fit_sphe.C)
    Double_t sphe_charge2  = 3425.95;
    Double_t sphe_std      = 500; // std dev of the first peak (not needed of deltaminus != 0)
    Double_t spe_max_val_at_time_cut = 20; // after `time_cut`, the signal cannot be higher than this
                                           // this allows to remove after pulses
    Double_t time_cut = 2000; // in ns seconds

    // coeficients to surround gaussian of 1 spe.
    // Gain             = (sphe_charge2 - sphe_charge)
    // spe's get events where charge < Gain*deltaplus  and charge < Gain/deltaminus
    // If deltaminus is set to zero, sphe_std*deltaplus will be used instead
    // This value can be checked with fit_sphe.C
    Double_t deltaplus  = 1.3;
    Double_t deltaminus = 1.5;


    // Not so common to change
    Double_t social_distance = 2; // demands that there is a minimum distance of social_distance * timeHigh between 2 consecutive peaks found
    string   method          = "derivative"; // `derivative` uses the 1th derivative of the waveform and a fixed threshold
                                             // `dynamic` evaluation of the baseline, search over moving average waveform
                                             // `fix` will not evaluate baseline and use raw threshold over moving average
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

    ADC_DATA sample; // will store the waveforms found

    string myname; // name of the declared class
    Bool_t derivate = false; // control if user has set to derivate
    Int_t kch; // current channel used
    Int_t nshow_finish = 100; // will show `nshow` debugging waveforms up
    Int_t nshow_start = 0; // starting from here
    Bool_t getstaticbase = false; // control if `static` method was set
    Double_t threshold = 0;

    Bool_t get_this_wvf = true;
    Bool_t get_this_charge = true;
    Int_t npts_wvf = 1000;
    Double_t *timeg = nullptr;

    string rootfile    = "analyzed.root";
    Int_t n_points;


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
    vector<Double_t> derivateMax; // Stores the max of the derivative
    vector<Double_t> derivateMaxFound; // Stores the max of the derivative before applying threshold

    vector<Double_t> smooted_wvf; // not needed to reset this
    vector<Double_t> denoise_wvf;

    vector<vector<Double_t>> selected_peaks; // these are for plotting
    vector<vector<Double_t>> selected_time;
    vector<Double_t> temp_selected_peaks;
    vector<Double_t> temp_selected_time;
    vector<Double_t> selected_charge;
    vector<Double_t> selected_charge_time;


    vector<Double_t> discardedTime; //store position of thed
    vector<Double_t> discardedPeak; //store position of thed
    vector<Double_t> discardedCharge; //store position of thed
    vector<Int_t> discarded_idx; //store idx (0, 1 or 2) of the peaks
    const char *selec_types[3] = {"Social distance", "Negative hits", "Too big"};
    TH1D *hdiscard = nullptr;
    vector<Double_t> mean_waveform;
    Int_t naverages = 0;
    // ____________________ ________________________________ ____________________ //

    vector<string> discard_reason = {"Social distance", "Negative hits", "Too big"};

    SPHE2(string m_name) : myname{m_name}{
      hbase        = new TH1D(Form("hbase_sphe_%s",myname.c_str()),"histogram for baseline",5*800,-400,400);
      hbase_smooth = new TH1D(Form("hbase_smooth_sphe_%s",myname.c_str()),"histogram for baseline smoothed",5*800,-400,400);
      hcharge      = new TH1D(Form("hcharge_sphe_%s",myname.c_str()),Form("hcharge_sphe_%s",myname.c_str()),50000,0,0);
      hdiscard     = new TH1D(Form("hdiscard_sphe_%s",myname.c_str()),Form("hdiscard_sphe_%s",myname.c_str()),3,0,3);



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
     * Clean the selection
     * Integrate peaks and store the waveform
     * Plot graphs for debugging
     * Save histogram
     */
    void giveMeSphe(){
      gROOT->SetBatch(kTRUE);

      hcharge->Reset();
      rootfile = filename + ".root";

      if(!z){ // load analyzer in case it is nullptr
        z = new ANALYZER(myname.c_str());
        z->dtime = dtime;
        z->setAnalyzer(rootfile);
      }
      if(!z->t1){ // otherwise, check if the ttree was loaded
        z->dtime = dtime;
        z->setAnalyzer(rootfile);
      }
      z->setChannel(Form("Ch%d.",channel));

      kch = z->kch;
      z->getWaveform(kch);
      n_points = z->n_points;
      smooted_wvf.resize(n_points);
      denoise_wvf.resize(n_points);
      fout = new TFile(Form("sphe_histograms_Ch%i.root",channel),"RECREATE");

      // ____________________ Setting up what needs to be set ____________________ //

      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
        twvf = new TTree("t1","mean waveforms");
      }
      else{ // in the case we are not taking waveform, I change this values for the same setup of integration
        mean_before = timeLow;
        mean_after = timeHigh;
      }
      // if the user set mean_before = 24 and mean_after = 100
      // we need to get  6 + 25 points + 1 = 32  (+1 because of the peak start)
      // if I found a point in around 100 ns, that is i = 25, I perform a while that goes
      // from i-mean_before = 19 up to i+mean_after = 50 (included) = 32 pts
      npts_wvf = (mean_before/dtime + mean_after/dtime) + 1;
      if (led_calibration) npts_wvf = (mean_after/dtime - mean_before/dtime);
      
      mean_waveform.clear();
      mean_waveform.resize(npts_wvf);
      naverages = 0;

      sample.Set_npts(npts_wvf);
      if (get_wave_form) twvf->Branch(Form("Ch%i",channel),&sample);


      nshow_start = nshow_range[0];
      nshow_finish = nshow_range[1];
      n_points = z->n_points;
      timeg = new Double_t[n_points];
      for(int i = 0; i < n_points; i++){
        timeg[i] = i*dtime;
      }
      if(led_calibration==true){
        method = "led";
      }
      else{
        // this might be too much, time will tell. Keep an eye in the RAM memory
        discardedTime.reserve(n_points);
        discardedPeak.reserve(n_points);
        discardedCharge.reserve(n_points);
        discarded_idx.reserve(n_points);
        peakPosition.reserve(n_points);
        peakMax.reserve(n_points);
        peaksFound.reserve(n_points);
        derivateMax.reserve(n_points);
        derivateMaxFound.reserve(n_points);
        temp_selected_peaks.reserve(n_points);
        temp_selected_time.reserve(n_points);
        selected_charge.reserve(n_points);
        selected_charge_time.reserve(n_points);
      }

      if(method == "static") getstaticbase = true;
      else if(method == "fix"){
      }
      else if(method == "derivative"){
        derivate = true;
        peaksRise.reserve(n_points);
        peaksCross.reserve(n_points);
      }

      Double_t delta = sphe_charge2 - sphe_charge;
      if(deltaminus == 0){
        delta1 = sphe_charge - deltaplus*sphe_std;
        delta2 = sphe_charge + deltaplus*sphe_std;
      }
      else{
        delta1 = delta/deltaminus;
        delta2 = delta*deltaplus;
      }

      // ____________________ ___________________________ ____________________ //



      Int_t nentries = z->nentries;
      if(just_a_test){nentries = just_this;}
      for(Int_t i = 0; i<nentries; i++){
        theGreatReset();
        z->getWaveform(i,kch); // get waveform by memory
        if(check_selection && z->ch[kch]->selection != 0) continue;
        
        if(i>nshow_finish && static_cast<Int_t>(i)%200==0){
          cout << i << " out of " << nentries << "\r" << flush;
        }
        processData();

        if(method == "static"){
          searchForPeaks();
        }
        else if(derivate){
          // Search for crossing points in the derivated waveform
          // positive crossing stored in peaksRise
          // negative crossing stored in peaksCross
          // from start to finish, discounting what was for getting integral with an extra point (not waveform, this filter is done later)
          Double_t start_of_search_offset = (timeLow>interactions*dtime) ? timeLow : interactions*dtime;
          Double_t finish_of_search_offset = (timeHigh>interactions*dtime) ? timeHigh : interactions*dtime;
          z->zeroCrossSearch(&smooted_wvf[0], peaksRise, peaksCross, start+start_of_search_offset+dtime, finish-finish_of_search_offset-dtime);
          derivateApplyThreshold();
        }
        cleanPeaks();
        integrateSignals(sample);
        if(snap()){
          drawMySamples();
          if(!led_calibration){
            cout << "Event " << z->ch[kch]->event << " total of peaks: " << peaksFound.size() << ", Valid = " << selected_charge.size() << endl;
            for(unsigned int j = 0; j < selected_charge.size(); j++){
              cout << "\t\t Charge = " << selected_charge[j] << " at " << selected_charge_time[j] << " ns " << endl;
            }
          }
        } // draw sample graphs
      }
      cout << "\n\n" << endl;


      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",kch));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
        Int_t ndiv = (naverages == 0 ? 1 : naverages);
        for(Int_t i = 0; i<npts_wvf; i++){
          mean_waveform[i] = mean_waveform[i]/ndiv;
        }

        TGraph *gmean = new TGraph(npts_wvf,timeg,&mean_waveform[0]);
        fwvf->WriteObject(gmean,"mean","TObject::kOverwrite");
        cout << "A total of " << naverages << " waveforms where found "<< endl;
      }

      fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),z->getIdx()),"TObject::kOverwrite");
      fout->WriteObject(hdiscard,"hdiscard","TObject::kOverwrite");


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
          temp_selected_peaks.clear();
          temp_selected_time.clear();
          selected_charge.clear();
          selected_charge_time.clear();
          discardedTime.clear();
          discardedPeak.clear();
          discardedCharge.clear();
          discarded_idx.clear();
          peaksFound.clear();
          derivateMax.clear();
          derivateMaxFound.clear();
        }
      
        denoise_wvf.clear();
        denoise_wvf.resize(n_points);
        smooted_wvf.clear();
        smooted_wvf.resize(n_points);
    }

    void processData(){
      Double_t *vec = nullptr;

      z->applyDenoise(filter, z->ch[kch]->wvf, &denoise_wvf[0]); // denoise is stored at denoise_wvf. If no filter, they are equal

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
        if(derivate_raw) z->differenciate(diff_multiplier,z->ch[kch]->wvf,&smooted_wvf[0]); // multiply by 1e3 so we can see something :)
        else z->differenciate(diff_multiplier,&denoise_wvf[0],&smooted_wvf[0]); // multiply by 1e3 so we can see something :)
        vec = &smooted_wvf[0];
      }
      else{
        vec = z->ch[kch]->wvf;
      }
      for(Int_t i = 0; i < ninter; i++){
        if (i == interactions-1 && method == "static") getstaticbase = true;
        smooted_wvf = movingAverage(vec, interactions, getstaticbase);
        vec = &smooted_wvf[0];
      }

    }

    bool goodSocialDistance(Int_t id2, Int_t id1){
      if(abs(id2-id1) <= social_distance*timeHigh/dtime){
        // if(snap()) cout << id2 << " " << id1 << " " << abs(id2-id1) << " " <<  social_distance*timeHigh/dtime << endl;
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

        return false;
      }

    }

    void derivateApplyThreshold(){
      Int_t ntotal = (int)peaksCross.size();
      if((int)peaksRise.size() < ntotal){ // there should be only pairs
        cout << "@@@@@@@@@@@@@@@@@@@@ THIS SHOULD NOT HAPPEN @@@@@@@@@@@@@@@@@@@@" << endl;
        cout << "\n\n\n\n\n\n\n\n\n\n" << endl;
        cout << "@@@@@@@@@@@@@@@@@@@@ THIS SHOULD NOT HAPPEN @@@@@@@@@@@@@@@@@@@@" << endl;
        ntotal = (int)peaksRise.size();
      }

      for(Int_t i = 0; i < ntotal; i++){
        Double_t crossPositive = peaksRise[i];
        Double_t candidatePosition = peaksCross[i]; // peaks previously selected

        Double_t dmax = z->getMaximum(crossPositive*dtime, candidatePosition*dtime, &smooted_wvf[0]);
        // if(snap()) cout << crossPositive*dtime << " " << candidatePosition*dtime << " " << dmax  << "\n";
        if(dmax < tolerance) {
          continue;
        }

        derivateMaxFound.push_back(dmax);
        peaksFound.push_back(peaksCross[i]) ;
      }

    }
    void cleanPeaks(){
      Int_t ntotal = (int)peaksFound.size();
      Bool_t wait_now = false;
      Double_t refWait = 0;// reference time that started waiting from a big pulse

      Bool_t discard_by_distance = false;
      Bool_t current_good_distance = false;

      for(Int_t i = 0; i < ntotal; i++){
        Double_t candidatePosition = peaksFound[i]; // peaks previously selected

        if(checkTooBig(wait_now, refWait, candidatePosition)){
          discard_by_distance = false;
          // if(snap()){
          //   // cout << wait_now << " " << refWait << " " << candidatePosition << " " << waiting << endl;
          // }
          fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 2);
          continue;
        }

        if(wait_now){ // I am waiting for a big pulse to finish
          if(candidatePosition*dtime >= (refWait+waiting)){
            // if(snap()){
            //   cout << "Got out of big pulse: " << candidatePosition << endl;
            // }
            wait_now = false;
          }
          else{
            fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 2);
            continue;
          }
        }

        //request social distance
        if( (ntotal > 1  && !(current_good_distance = goodSocialDistance(peaksFound[i+1],candidatePosition))) || discard_by_distance){
          // if(snap()) cout << discard_by_distance << endl;
          if (!current_good_distance) discard_by_distance = true;
          else discard_by_distance = false;

          // discarding current peak by social
          fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 0);
          continue;
        }
        peakPosition.push_back(peaksFound[i]);
        derivateMax.push_back(derivateMaxFound[i]);
      }
    }

    void fill_discarded(Int_t pidx, Double_t val, Double_t charge, Int_t discardidx){
      if(snap())
      {
        discardedTime.push_back(pidx*dtime);
        discardedPeak.push_back(val);
        discardedCharge.push_back(charge);
        discarded_idx.push_back(discardidx);
      }
      hdiscard->Fill(selec_types[discardidx], 1);
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
      Int_t n = n_points;

      vector<Double_t> res(n,0);
      if(myinte==0) { // nothing to do
        for (Int_t i = 0; i < n; i++) {
          res[i] = v[i];
        }
        return res;
      }
      else{
        interactions = myinte;
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

    void integrateSignals(ADC_DATA &sample){
      unsigned int npeaks = peakPosition.size();
      if (led_calibration) npeaks = 1;
      for(unsigned int i = 0; i < npeaks; i++){
        Int_t peakPosIdx = peakPosition[i]; // position of the peak as idx ( don't worry for led )
        temp_selected_peaks.clear();
        temp_selected_time.clear();
        getIntegral(peakPosIdx, sample); // If get_mean_wvf is set to false, there is no problem!
        Double_t charge = z->temp_charge;
        Double_t peak = z->temp_max;
        if(get_this_charge){
          hcharge->Fill(charge);
          if(snap())
          {
            selected_charge.push_back(charge);
            selected_charge_time.push_back(peakPosIdx*dtime);
            selected_peaks.push_back(temp_selected_peaks);
            selected_time.push_back(temp_selected_time);
          }
        }
        else{
          continue; // If I didn't get the charge, I will not take the waveform
        }
        if(get_this_wvf){
          sample.selection = 0;
          if(charge>=delta1 && charge <= delta2){
            sample.selection = 1;
            for(Int_t j = 0; j < npts_wvf; j++){
              mean_waveform[j] += sample.wvf[j];
              if(j*dtime >= time_cut && sample.wvf[j] >= spe_max_val_at_time_cut){
                sample.selection = 2;
                for(Int_t k = j; k >= 0; k--){
                  mean_waveform[k] -= sample.wvf[k];
                }
                break;
              }
            }
            if(sample.selection == 1) naverages+=1;
          }
          sample.peak = peak;
          sample.peakpos = peakPosIdx*dtime;
          sample.charge = charge;
          sample.event = z->currentEvent;
          sample.fprompt = derivateMax[i];
          twvf->Fill();
        }
      }
    }

    void getIntegral(Int_t peakPosIdx, ADC_DATA &sample){
      get_this_charge = true;
      if(get_wave_form) get_this_wvf = true; // in the case not, I dont take the waveform
      else get_this_wvf = false;
      Double_t from = peakPosIdx - mean_before/dtime;
      Double_t to = peakPosIdx + mean_after/dtime;
      
      if(from < 0 || to >= n_points){ // check if boundaries are being respected
        from = peakPosIdx - timeLow/dtime;
        to = peakPosIdx + timeHigh/dtime;
        if(get_wave_form) get_this_wvf = false; // in the case not, I dont take the waveform
      }
      Int_t iwvf = 0;

      Double_t res = 0;
      Double_t max = -1e12;
      Int_t negativeHits = 0;
      Double_t integralfrom = peakPosIdx - timeLow/dtime;
      Double_t integralto = peakPosIdx + timeHigh/dtime;

      if(led_calibration){
        from = 0;
        to = n_points-1;
        integralfrom = start/dtime;
        integralto = finish/dtime;
        get_this_wvf = true;
      }
      for(Int_t i = from; i <= to; i++, iwvf++){
        Double_t val = denoise_wvf[i];
        if(get_this_wvf){
          sample.wvf[iwvf] = z->ch[kch]->wvf[i];
        }
        if(i >= integralfrom && i <= integralto){
          if (val < lowerThreshold){
            negativeHits++;
          }
          
          if(negativeHits >= maxHits){
            get_this_wvf = false;
            get_this_charge = false;
            fill_discarded(peakPosIdx, smooted_wvf[peakPosIdx], res*dtime, 1);
            break;
          }
          if(snap())
          {
            temp_selected_peaks.push_back(val);
            temp_selected_time.push_back(i*dtime);
          }
          if(withfilter == false){ val = z->ch[kch]->wvf[i]; }
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
      if(z->currentEvent >= nshow_start && z->currentEvent < nshow_finish && !led_calibration){
        return true;
      }
      else{
        return false;
      }
      return false; // just to be sure oO
    }
    void drawMySamples(){
      unsigned int npeaksselected = selected_charge.size();
      unsigned int npeakstotal = peaksFound.size();
      string sampleName = "ev_" + to_string(z->currentEvent)+"_"+to_string(npeaksselected)+"_of_"+to_string(npeakstotal);

      TGraph g_smooth(n_points, timeg, &smooted_wvf[0]);
      TGraph g_normal(n_points, timeg, &denoise_wvf[0]);
      // TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str(),1920,0,700,500);
      // this is not working when saving
      TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str());


      c1->cd(1);
      g_smooth.SetLineColor(kRed);
      g_smooth.SetLineWidth(3);

      g_normal.SetLineColor(kBlue);
      g_normal.SetTitle(" ");

      g_normal.SetEditable(kFALSE);
      g_smooth.SetEditable(kFALSE);

      g_normal.Draw("AL");
      g_smooth.Draw("L SAME");

      if (derivate || getstaticbase ) {
        threshold = tolerance;
        mean = 0;
      }

      TLine lmean(start,mean,finish,mean);
      TLine ldev(start,mean+threshold,finish,mean+threshold);

      lmean.SetLineColor(kGreen);
      ldev.SetLineColor(kGreen);

      lmean.SetLineWidth(2);
      ldev.SetLineWidth(2);

      lmean.Draw();
      ldev.Draw();
      Int_t n = selected_charge.size();
      vector<TGraph*> g_selected(n);
      for(Int_t i = 0; i < n; i++){
        g_selected[i] = new TGraph(selected_time[i].size(), &selected_time[i][0], &selected_peaks[i][0]);
        g_selected[i]->SetMarkerColor(kBlack);
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << selected_charge[i];
        string gtitle = "charge = " + stream.str();
        g_selected[i]->SetTitle(gtitle.c_str());
        g_selected[i]->SetEditable(kFALSE);
        g_selected[i]->Draw("P* SAME");
      }

      Int_t ndis = discardedPeak.size();
      TGraph *g_discarded = nullptr;
      Int_t marker_style = 21;
      Color_t marker_color = kRed;
      for(Int_t i = 0; i < ndis; i++){
        g_discarded = new TGraph(1, &discardedTime[i], &discardedPeak[i]);
        stringstream stream_discarded;
        stream_discarded << std::fixed << std::setprecision(2) << discardedCharge[i];
        string gtitle_discarded = "charge = " + stream_discarded.str();
        g_discarded->SetTitle(gtitle_discarded.c_str());
        if(discarded_idx[i] == 0){ // social distance
          marker_style = 22; // triangle
          marker_color = kYellow;
        }
        else if(discarded_idx[i] == 1){ // negative hits
          marker_style = 29; // filled star
          marker_color = kGreen;
        }
        else if(discarded_idx[i] == 2){ // Too big
          marker_style = 21; // filled square
          marker_color = kRed;
        }
        g_discarded->SetMarkerSize(2);
        g_discarded->SetMarkerStyle(marker_style);
        g_discarded->SetMarkerColor(marker_color);
        g_discarded->SetEditable(kFALSE);
        g_discarded->Draw("SAME P");
      }

      // if(derivate) z->drawZeroCrossingLines(peaksCross, peaksRise, c1,0,tolerance);
      fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");
      for(Int_t i = 0; i < n; i++){
        delete g_selected[i];
      }
      g_selected.clear();

      delete c1;


    }



};



























class MeanSignal{
  
  public:
  
    Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits
    Int_t n_points = memorydepth;
  
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
      vector<OLD_ADC_DATA> ch(channels.size());
      vector<TBranch*> bch(channels.size());
      for(Int_t k = 0; k<(int)channels.size();k++){
        bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
        bch[k]->SetAddress(&ch[k]);
      }
      Int_t nentries = t1->GetEntries();
      vector<Double_t> norm(channels.size(),0);
  
      vector<vector<Double_t>> avg(channels.size(),std::vector<Double_t>(n_points,0));
      vector<vector<Double_t>> avgn(channels.size(),std::vector<Double_t>(n_points,0));
      vector<Double_t> time(n_points);
      for(Int_t j = 0; j<n_points; j++){
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
              for(Int_t j = 0; j<n_points; j++){
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
        for(Int_t j = 0; j<n_points; j++){
          avg[k][j]=avg[k][j]/norm[k];
          avgn[k][j]=avg[k][j]/maxvalue[k];
        }
        cout << "Ch" << channels[k] << " total waveforms = " << norm[k] << endl;
        gavg[k] = new TGraph(n_points,&time[0],&avg[k][0]);
        gavgn[k] = new TGraph(n_points,&time[0],&avgn[k][0]);
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
    Int_t n_points = memorydepth;

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
      vector<OLD_ADC_DATA> ch(channels.size());
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
      for(Int_t i = 0; i<n_points; i++){
        if(i>=intmin/dtime && i<(intmin+integrations[index_ints])/dtime){
          refCharge+=wvf[i]*dtime;
          if(i==n_points-1){
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
    Int_t n_points = memorydepth;
  
  
    vector<double> delay_line(Double_t v[], Int_t delay_time){
      vector<double> res(n_points);
      for(int i=0; i<n_points; i++){
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
      vector<vector<Double_t>> peaks(channels.size(),std::vector<Double_t>(n_points,0));
      vector<vector<Double_t>> shifted(channels.size(),std::vector<Double_t>(n_points,0));
      vector<TGraph *> gsample(nshow*nchannels);
      vector<TGraph *> gshifted(nshow*nchannels);
      vector<TH1D*> htd(nchannels);
      vector<Int_t> posaux(nchannels);
      vector<vector<Double_t>> position(nchannels);
  
      string rootfile = filename + ".root";
      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");
  
      TFile *fout = new TFile("time_distribuction.root","RECREATE");
  
      vector<OLD_ADC_DATA>  ch(nchannels);
      vector<TBranch *> bch(nchannels);
      for(Int_t k = 0; k<nchannels;k++){
        bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
        bch[k]->SetAddress(&ch[k]);
        htd[k] = new TH1D(Form("h_ch%i",channels[k]),Form("h_ch%i",channels[k]),n_points,0,dtime*n_points);
      }
      for(Int_t j = 0; j<n_points; j++){
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
          for(Int_t j = 1; j<n_points-1; j++){
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
            gsample[nchannels*i+k] = new TGraph(n_points,&gtime[0],&ch[k].wvf[0]);
            gshifted[nchannels*i+k] = new TGraph(n_points,&gtime[0],&shifted[k][0]);
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
        if(i<n_points){
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
  
      for(Int_t i = 0; i<n_points; i++){
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
              
