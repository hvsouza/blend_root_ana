class OUTPUT{
public:
  Double_t charge;
  Double_t slow;
    
  OUTPUT();
  const char *tobranch = Form("charge/D:slow/D");
};

OUTPUT::OUTPUT(){
  charge = 0;
  slow = 0;
}

Double_t sphe_charge = 0.3; 

const Int_t nparameters = 1;

vector<Double_t> volts = {0,10};
Int_t n = volts.size();
Int_t n_bins = 0;
TFile *f;
TTree *t1;
  
vector<OUTPUT> out(n);
vector<TBranch*> br(n);
  
vector<TH1D*> hcharge(2);
vector<TH1D*> hslow(2);
TH1D* hfit;
TH1D* hfitCharge;
 
vector<Double_t> y;
vector<Double_t> z;

Int_t thisone = 0;
string type = "charge";
// string type = "slow";

Double_t fit_start = 0; // in charge
Double_t fit_finish = 5300; // 1e12 if no maximum
Double_t fit_start_bin = 0;
Double_t fit_finish_bin = 1e12;

vector<Double_t> scaling_factor;
vector<Double_t> fchisq;

vector<Double_t> tmp_scaling_factor;
vector<Double_t> tmp_fchisq;

void makeSpectrum(Double_t &val, Double_t par[]){
  Double_t res;
  for(Int_t j = 0; j<t1->GetEntries(); j++){
    br[0]->GetEvent(j);
    res = val*par[0];
    hfit->Fill(res/sphe_charge);
  }
}
Int_t manytimes = 0;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chisquare
    Double_t chisq = 0;
    Double_t delta =0;
  
  // makeSpectrum_Noise(data[0].detNQ,param);
    hfit->Reset();

    if(type == "charge") makeSpectrum(out[0].charge,par);
    else makeSpectrum(out[0].slow,par);

    for(Int_t i=0; i<n_bins; i++) {
      y[i] = (hfit->GetBinContent(i+1));
      delta  = (y[i]-z[i]);
      // cout << y[i] << " " << z[i] <<  " " << i << " " << fit_finish_bin << endl;
      if(z[i]>0 && i>fit_start_bin && i < fit_finish_bin)chisq += delta*delta/z[i];
    }
    f = chisq;
    fchisq.push_back(f);
    scaling_factor.push_back(par[0]);
    manytimes++;
    cout << manytimes << ", " << f;
    
    for(Int_t i = 0; i < nparameters; i++){
      cout << ", " << par[i];
    }
    cout << "\n";
}


void findMin(vector<Double_t> x, Int_t &ref){
  Int_t n = x.size();
  Int_t minpos = 0;
  Double_t min = 1e12;
  for(Int_t i = ref; i<n; i++){
    if(x[i]<=min){
      min = x[i];
      minpos = i;
    }
  }
  // cout << minpos << " " << min << " ";
  ref = minpos;
  // return minpos;
}
void theGreatReset(vector<Double_t> &x, vector<Double_t> &y){
  Int_t n = y.size();
  Int_t minpos = 0;
  Int_t ref = 0;
  Double_t temp = 0;
  for(Int_t i = 0; i<n; i++){
    findMin(x,ref);
    temp = x[ref];
    x[ref] = x[i];
    x[i] = temp;
    temp = y[ref];
    y[ref] = y[i];
    y[i] = temp;
    ref = i+1;
    // cout << x[i] << " " << y[i] << endl;
  }
}


void fit_C2_histograms(){
  Int_t n_bins_charge = 150;
  Int_t n_bins_slow = 50;
  Double_t hist_max = 0;
  Double_t hist_max_charge = 26700;
  Double_t hist_max_slow = 3300;

  if(type=="charge"){
    n_bins = n_bins_charge;
    hist_max = hist_max_charge;
  }
  else{
    n_bins = n_bins_slow;
    hist_max = hist_max_slow;
  }
  y.resize(n_bins);
  z.resize(n_bins);
 


  hfit = new TH1D("hfit","hfit",n_bins,0,hist_max);
  hfitCharge = new TH1D("hfitCharge","hfitCharge",n_bins,0,hist_max);
  for(Int_t k = 1; k<=n_bins; k++){
    if(hfit->GetBinCenter(k)>fit_start){
      fit_start_bin = k;
      break;
    }
  }
  for(Int_t k = 1; k<=n_bins; k++){
    if(hfit->GetBinCenter(k)>fit_finish){
      fit_finish_bin = k;
      break;
    }
  }
  // cout << fit_start << " at bin " << fit_start_bin << " corresponde to " << hfit->GetBinCenter(fit_start_bin) << endl;
  
  f = new TFile("C2_output_charges.root","READ");
  t1 = (TTree*)f->Get("t1");
  for (Int_t i = 0; i<n; i++){
    br[i] = t1->GetBranch(Form("E%.0f",volts[i]));
    br[i]->SetAddress(&out[i]);
    if(type=="charge") hcharge[i] = new TH1D(Form("hcharge_E%.0f",volts[i]),Form("hcharge_E%.0f",volts[i]),n_bins_charge,0,hist_max_charge);
    else hcharge[i] = new TH1D(Form("hcharge_E%.0f",volts[i]),Form("hcharge_E%.0f",volts[i]),n_bins_slow,0,hist_max_slow);
    // hslow[i] = new TH1D(Form("hslow_E%.0f",volts[i]),Form("hslow_E%.0f",volts[i]),n_bins_slow,0,hist_max_slow);
    for(Int_t j = 0; j<t1->GetEntries(); j++){
      br[i]->GetEvent(j);
      if(type=="charge") hcharge[i]->Fill(out[i].charge/sphe_charge);
      else hcharge[i]->Fill(out[i].slow/sphe_charge);
    }
  }

  for(Int_t k = 0; k<n_bins; k++){
    z[k] = hcharge[1]->GetBinContent(k+1);
  }
  // hcharge[0]->Draw("");
  // hcharge[1]->Draw("same");

  
  


  TMinuit *gMinuit_charge = new TMinuit(nparameters);  //initialize TMinuit with a maximum of 3 params
  gMinuit_charge->SetFCN(fcn); // setting fcn 
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  gMinuit_charge->mnexcm("SET ERR", arglist ,1,ierflg);

  Double_t par1 = 0.8;
  Double_t par2 = 0.75;
  // Set starting values and step sizes for parameters
  Double_t vstart[nparameters] = {par1};
  Double_t step[nparameters] = {abs(par1)/10.};
  gMinuit_charge->mnparm(0, "scale", vstart[0], step[0], 0,2,ierflg);
  // gMinuit_charge->mnparm(1, "norm", vstart[1], step[1], 0.1,20,ierflg);

  // Now ready for minimization step
  arglist[0] = 0;
  arglist[1] = 1000;
  arglist[2] = 0.2;
  arglist[3] = 1;
  

  gMinuit->mnexcm("SCAn", arglist ,4,ierflg);

  // Now ready for minimization step
  arglist[0] = 10000000;
  arglist[1] = 1000;

 
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  
  // gMinuit->mnexcm("IMPRO", arglist ,1,ierflg);

  
  Double_t param[nparameters];
  Double_t Erparam[nparameters];
  
  for(Int_t i = 0; i < nparameters; i++){
    gMinuit->GetParameter(i,param[i],Erparam[i]);
    cout << param[i] << " +/- " << Erparam[i] << endl;
  }
  
  hfit->Reset();


  if(type=="charge") makeSpectrum(out[0].charge,param);
  else makeSpectrum(out[0].slow,param);
  

  TCanvas *c1 = new TCanvas();
  hfitCharge = (TH1D*)hfit->Clone("hfitCharge");
  hfitCharge->SetLineWidth(2);
  hfitCharge->SetLineColor(kBlack);

  hfitCharge->SetNameTitle("HV Off reweighted","HV Off reweighted");
  hcharge[0]->SetNameTitle("HV Off","HV Off");
  hcharge[1]->SetNameTitle("HV = 10 kV","HV = 10 kV");


  hcharge[0]->GetXaxis()->SetTitle("Photo-electrons");
  // hcharge[0]->GetXaxis()->SetTitle("Photo-electrons from slow comp.");
  hcharge[0]->GetYaxis()->SetTitle("# of events");
  // hs->Add(hfitNQ,"hist sames");
  // hs->Draw("nostack");
  hcharge[0]->Draw("hist");
  hcharge[1]->Draw("hist same");
  hcharge[1]->SetLineWidth(2);
  hcharge[1]->SetLineColor(kRed);
  hfitCharge->Draw("hist same");
  hcharge[0]->SetLineWidth(2);
  hcharge[0]->SetLineColor(kBlue);
  gPad->SetLogy(1);

  c1->BuildLegend();
  
  hfit->Reset();


  theGreatReset(scaling_factor,fchisq);
  TCanvas *c2 = new TCanvas();
  TGraph *gchiq = new TGraph(fchisq.size(),&scaling_factor[0],&fchisq[0]);
  gchiq->Draw("APL");
  gchiq->SetMarkerStyle(21);
  gchiq->SetMarkerSize(0.5);
  
  gchiq->GetXaxis()->SetRangeUser(0.18,1.02);
  gchiq->GetXaxis()->SetTitle("Scaling factor");
  gchiq->GetYaxis()->SetTitle("#chi^{2}");
  




  
}
