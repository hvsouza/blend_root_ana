#define memorydepth 5000
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

Double_t mysignals(Double_t *x, Double_t *par){
  Double_t xx = x[0];
  if(xx>par[0]){
    return par[2]*exp(-(xx-par[0])/par[1]);
  }
  else{
    return 0;
  }
}

TH1D *hbase_smooth = new TH1D("hbase_smooth","hbase_smooth",400,0,0);
void smoothWithMovingAvarage(vector<Double_t> &peak_smooth,vector<Double_t> &v,Int_t interactions){

  Int_t midpoint = 0;
  Int_t width = 0;
  Double_t sum = 0;
  Int_t n = v.size();
  if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
    midpoint = interactions/2+1;    //midpoint will be 5 here
    width = interactions+1;
  }
  else{
    midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
    width = interactions;
  }
    
    
  for(Int_t i = 0; i < n; i++){
    if(interactions==0){
      peak_smooth.push_back(v.at(i));
      hbase_smooth->Fill(peak_smooth.at(i));
      continue;
    }
    if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
      peak_smooth.push_back(v.at(i));
    }
    else{
      for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
        sum = sum+v.at(j);
        //                 cout << sum << endl;
      }
      peak_smooth.push_back(sum/width);
      hbase_smooth->Fill(peak_smooth.at(i));
          
      
    }
    
    
        
    sum=0;
  }
    
    
}


vector<double> delay_line(vector<double> v, double delay_time){
    if(delay_time==0) return v;
    vector<double> res(v.size());
    for(int i=0; i<v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
    }
    return res;
}

void testing_offset(){

  vector<Double_t> raw(memorydepth);
  vector<Double_t> signal(memorydepth);
  vector<Double_t> fder(memorydepth);
  vector<Double_t> fder_ma;
  vector<Double_t> t(memorydepth);
  Double_t offset = 2431;
  TH1D *hbase = new TH1D("hbase","hbase",1000,0,1000*4);

  TFile *f1 = new TFile("/home/henrique/Documents/ADC_data/test_SuperCell/20211218/run8_31V50_400ADC_Ch0/analyzed.root","READ");
  TTree *t1 = (TTree*)f1->Get("t1");
  vector<Int_t> channels = {1};
  vector<ADC_DATA> ch(channels.size());
  vector<TBranch*> bch(channels.size());
  for(Int_t k = 0; k<channels.size();k++){
    bch[k] = t1->GetBranch(Form("Ch%i",channels[k]));
    bch[k]->SetAddress(&ch[k]);
  }
  

  // creating fake expo like signalsignal
  gRandom->SetSeed(0);
  Int_t n_peaks = 5;
  Double_t starting_time = -200;
  Double_t ending_time = 8000;
  vector<TF1*> f(n_peaks);
  vector<Double_t> pktime(n_peaks);
  vector<Double_t> pkamp(n_peaks);

  
  
  for(Int_t i = 0; i<n_peaks; i++){
    f[i] = new TF1("f",mysignals,0,memorydepth*4,3);
    pktime[i] = gRandom->Uniform(starting_time,ending_time);
    pkamp[i] = gRandom->Uniform(10,12);
    f[i]->SetParameters(pktime[i],200,pkamp[i]);
  }

  vector<Double_t> found_peaks;
  bch[0]->GetEvent(14);
  for(Int_t i = 0; i<memorydepth; i++){
    for(Int_t j = 0; j<n_peaks; j++){
      raw[i] += f[j]->Eval(i*4);
    }
    raw[i] = raw[i]+gRandom->Gaus(0,2)+offset;
    signal[i] = ch[0].wvf[i];
    t[i] = i*4;
  }

  // DENOISE dn;
  // dn.TV1D_denoise<Double_t>(&raw[0],&signal[0],memorydepth,3);

  for(Int_t i = 0; i<2500; i++){
    hbase->SetBinContent(i+1,signal[i]);
  }


  WIENER wbase("wbase",4,250,1e-9,1e6,2500);
  wbase.fft(hbase);
  

  // TGraph *graw = new TGraph(memorydepth,&t[0],&raw[0]);
  TGraph *g = new TGraph(memorydepth,&t[0],&signal[0]);
  // graw->SetLineColor(kGreen);
  // graw->Draw("AL");
  g->Draw("AL SAME");

  TCanvas *c2 = new TCanvas();
  wbase.hfft->Draw();

  cout << wbase.hfft->GetBinContent(1) << endl;
  for(Int_t i = 0; i<memorydepth; i++){
    signal[i] = signal[i] - wbase.hfft->GetBinContent(1);
    // if(i>0){
      // fder[i] = signal[i]-signal[i-1];
    // }
  }
  fder = delay_line(signal,46);


  smoothWithMovingAvarage(fder_ma,fder,55);

  Double_t rms = 0;
  Double_t max = 1e-12;
  Double_t nrms = 0;
  for(Int_t i = 3800/4; i<4200/4; i++){
    if(fder[i]>max) max = fder_ma[i];
  }
  Double_t maxnoise = 1e-12;  
  for(Int_t i = 6000/4; i<14000/4; i++){
    rms += fder_ma[i]*fder_ma[i];
    nrms+=1;
   if(fder_ma[i]>maxnoise) maxnoise = fder_ma[i]; 
  }
  rms = sqrt(rms/nrms);

  cout << "SNR = " << max << "/" << rms << " = " << max/rms << endl;
  cout << "SNR = " << max << "/" << maxnoise << " = " << max/maxnoise << endl;

  TCanvas *c3 = new TCanvas();
  TGraph *g2 = new TGraph(memorydepth,&t[0],&fder[0]);
  TGraph *g2ma = new TGraph(memorydepth,&t[0],&fder_ma[0]);
  g2->SetLineColor(kGreen);
  g2->Draw("A L");
  g2ma->Draw("L SAME");

  


  TCanvas *c4 = new TCanvas();
  hbase_smooth->Draw();


  
}
