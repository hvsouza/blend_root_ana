#define memorydepth 36864
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void read(vector<Double_t> &time, vector<Double_t> &v, string filename){

  ifstream f;
  f.open(filename.c_str(),ios::in);
  if(f.good() && f.is_open()){ // Ok
    cout << "Reading file " << filename << " ... " << endl;
  }
  else{ 
    cout << "File " << filename << " did not open!!" << endl;
    return;
  }
  // reading headers
  string head;
  getline(f,head);
  getline(f,head);
  Double_t t, s;
  Int_t aux = 0;
  while(!f.fail()){
    f >> t >> s;
    if(f.bad() || f.fail()) break;
    time[aux] = t;
    v[aux] = s;
    aux++;
    if(aux>memorydepth){
      cout << "@@@@@@@@@@@@@@@@" << endl;
    }
  }
}

void denoise(){
  vector<Double_t> time(memorydepth);
  vector<Double_t> source(memorydepth);
  vector<Double_t> noise(memorydepth);

  read(time,source,"source_mono.dat");
  read(time,noise,"bg_mono.dat");

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,3);
  c1->cd(1);
  TGraph *gsource = new TGraph(memorydepth,&time[0],&source[0]);
  gsource->SetLineColor(kBlue);
  gsource->Draw("AL");
  c1->cd(2);
  TGraph *gnoise = new TGraph(memorydepth,&time[0],&noise[0]);
  gnoise->SetLineColor(kRed);
  gnoise->Draw("AL");


  WIENER wf;
  Int_t filter_size = 2000;
  vector<Double_t> w(filter_size);
  vector<Double_t> res(memorydepth);
  w = wf.create_filter(source,noise,filter_size);
  res = wf.filter(source,w);
  
  c1->cd(3);
  TGraph *gdenoise = new TGraph(memorydepth,&time[0],&res[0]);
  gdenoise->SetLineColor(kViolet-1);
  gdenoise->Draw("AL");

  ofstream fout;
  fout.open("source_filter.dat",ios::out);

  fout << "; Sample Rate 48000 \n";
  fout << "; Channels 1 \n";
  for(Int_t i =0; i<memorydepth; i++){
    fout << time[i] << " " << res[i] << endl;
  }
  fout.close();


 
  
}
