
class DATA{
public:
  Int_t memorydepth;
  vector<Double_t> x;
  vector<Double_t> y;
};
vector<string> infos = {"all_on","SiPM_off","Analog_off","Fiber_disconnected","Coherent_off","Cable_off_coherent","Cable_off_scope"};

vector<Color_t> mycolors = {kBlue,kGreen+2,kBlack,kYellow+3,kRed,kMagenta,kCyan+1}; 

void setCommon(TGraph *g){
  g->GetXaxis()->SetTitleOffset(1.1);
  // g->SetLineColor(kBlue);
  g->SetLineWidth(2);
}
void setWvfGraph(TGraph *g){
  g->GetXaxis()->SetTitle("Time (#mus)");
  g->GetYaxis()->SetTitle("Amplitude (V)");
  setCommon(g);
}
void setFFTGraph(TGraph *g){
  g->GetXaxis()->SetTitle("Frequency (kHz)");
  g->GetYaxis()->SetTitle("Amplitude (dB)");
  setCommon(g);
}



void readData(ifstream &f,DATA &data, Bool_t notFFT = true){

  string val;
  Double_t xtemp;
  Double_t ytemp;
  Double_t xref=0;
  Double_t timeScale = 1e6;
  // LECROYWaveRunner�R,21153,Waveform
  // Segments,1,SegmentSize,100002
  // Segment,TrigTime,TimeSinceSegment1
  // #1,16-Mar-2022 14:45:39,0                 
  // Time,Ampl
  getline(f,val);
  getline(f,val,',');
  getline(f,val,',');
  getline(f,val,',');
  getline(f,val);
  data.memorydepth = stoi(val);
  getline(f,val);
  getline(f,val);
  getline(f,val);

  // if(notFFT==false){
    // data.memorydepth=data.memorydepth-1;
  // }
   
  // data.x.resize(data.memorydepth);
  // data.y.resize(data.memorydepth);
  
  getline(f,val,',');
  xtemp = stod(val)*timeScale;
  xref = xtemp;
  // cout << xref << endl;
  getline(f,val);
  ytemp = stod(val);
  
  // if(notFFT){
    data.x.push_back(xtemp-xref);
    data.y.push_back(ytemp);
  // }   
  while(!f.fail()){
    for(Int_t i = 1; i<data.memorydepth; i++){
      getline(f,val,',');
      if(f.bad() || f.fail()){
        break;
      }
      xtemp = stod(val);
      if(notFFT) xtemp = xtemp*timeScale;
      else xtemp = xtemp*1e-3;
      if(i==0) xref = xtemp;
      getline(f,val);
      ytemp = stod(val);
      
      data.x.push_back(xtemp-xref);
      data.y.push_back(ytemp);
      // cout << data.y[i] << " " << data.x[i] << endl;
    }
    if(f.bad() || f.fail()){
      break;
    }
  }

}
void plot_graphs(){
  Int_t n = 7;
  vector<DATA> data(n);
  vector<DATA> data2(n);
  vector<ifstream> fwvf(n);
  vector<ifstream> ffft(n);
  vector<TGraph*> gwvf(n);
  vector<TGraph*> gfft(n);

  vector<TCanvas *> cwvf(n);
  vector<TCanvas *> cfft(n);
  
  for(Int_t i = 0; i<n; i++){
    fwvf[i].open(Form("%s/C4--100V-2000V--00000.txt",infos[i].c_str()),ios::in);
    ffft[i].open(Form("%s/SpecAn--100V-2000V--00000.txt",infos[i].c_str()),ios::in);
    if(fwvf[i].good() && fwvf[i].is_open()){ // Ok
    }
    else{ 
      cout << "File did not open!!" << endl;
      return;      
    }
    readData(fwvf[i],data[i]);    
    readData(ffft[i],data2[i],false);    

    gwvf[i] = new TGraph(data[i].memorydepth, &data[i].x[0],&data[i].y[0]);
    gfft[i] = new TGraph(data2[i].memorydepth, &data2[i].x[0],&data2[i].y[0]);
    gwvf[i]->SetNameTitle(Form("wvf - %s",infos[i].c_str()),Form("wvf - %s",infos[i].c_str()));
    gfft[i]->SetNameTitle(Form("fft - %s",infos[i].c_str()),Form("fft - %s",infos[i].c_str()));
    setWvfGraph(gwvf[i]);
    setFFTGraph(gfft[i]);
    gwvf[i]->SetLineColor(mycolors[i]);
    gfft[i]->SetLineColor(mycolors[i]);
    gwvf[i]->GetXaxis()->SetRangeUser(0,200000);
    gfft[i]->GetXaxis()->SetRangeUser(0,100);

    cwvf[i] = new TCanvas(Form("cwvf%d",i), Form("cwvf%d",i),1920,0,1920,1080);
    cfft[i] = new TCanvas(Form("cfft%d",i), Form("cfft%d",i),1920,0,1920,1080);
    cwvf[i]->cd();
    gwvf[i]->Draw("AL");
    cwvf[i]->BuildLegend();
    cfft[i]->cd();
    gfft[i]->Draw("AL");
    cfft[i]->BuildLegend();
    
  }

  TCanvas *call = new TCanvas("call","call",1920,0,1920,1080);
  TMultiGraph *gm = new TMultiGraph("gm","gm");
  gm->Add(gfft[0]);
  gm->Add(gfft[1]);
  gm->Add(gfft[3]);
  gm->Add(gfft[4]);
  gm->Add(gfft[5]);
  gm->Add(gfft[6]);
  gm->GetXaxis()->SetTitle("Frequency (kHz)");
  gm->GetYaxis()->SetTitle("Amplitude (dB)");
    
  gm->Draw("A L");
  call->BuildLegend();
}
