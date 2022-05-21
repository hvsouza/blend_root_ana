#define memorydepth 36864
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void read_with_a_twist(vector<Double_t> &time, vector<Double_t> &v, string filename, Int_t file_size = memorydepth){

  
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
  string head1;
  string head2;
  getline(f,head1);
  getline(f,head2);
  Double_t t, s;
  Int_t aux = 0;
  ofstream fout;
  fout.open("source_gaus.dat",ios::out);
  fout << head1 << endl;
  fout << head2 << endl;

  for(Int_t i = 0; i<file_size; i++){
    f >> t >> s;
    if(f.bad() || f.fail()) break;
    time[aux] = t;
    v[aux] = s+gRandom->Gaus(0,0.2);
    fout << t << " " << v[aux] << endl;
    aux++;
    if(aux>file_size){
      cout << "@@@@@@@@@@@@@@@@" << endl;
    }
  }

  }
  
void read(vector<Double_t> &time, vector<Double_t> &v, string filename, Int_t file_size = memorydepth){

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
  for(Int_t i = 0; i<file_size; i++){
    f >> t >> s;
    if(f.bad() || f.fail()) break;
    time[aux] = t;
    v[aux] = s;
    aux++;
    if(aux>file_size){
      cout << "@@@@@@@@@@@@@@@@" << endl;
    }
  }
}

void denoise_with_reference(){
  vector<Double_t> time(memorydepth);
  vector<Double_t> reference(memorydepth);
  vector<Double_t> signal(memorydepth);
  Int_t n2 = 204800;
  vector<Double_t> time2(n2);  
  vector<Double_t> speech(n2);  
  read(time,reference,"source_mono.dat");
  read_with_a_twist(time,signal,"source_mono.dat");
  read_with_a_twist(time2,speech,"speech2_mono.dat",n2);
  // read(time,signal,"bg_mono.dat");

  TCanvas *c1 = new TCanvas();
  c1->Divide(1,2);
  c1->cd(1);
  TGraph *gsource = new TGraph(memorydepth,&time[0],&reference[0]);
  gsource->SetLineColor(kBlue);
  gsource->Draw("AL");
  c1->cd(2);
  TGraph *gsignal = new TGraph(memorydepth,&time[0],&signal[0]);
  gsignal->SetLineColor(kRed);
  gsignal->Draw("AL");


  WIENER wf;
  Int_t filter_size = 200;
  vector<Double_t> w(filter_size);
  vector<Double_t> res(memorydepth);
  w = wf.create_filter_from_reference(signal,reference,filter_size);
  res = wf.filter(signal,w);
  vector<Double_t> res2(n2);
  res2 = wf.filter(speech,w);
  
  c1->cd(2);
  TGraph *gdenoise = new TGraph(memorydepth,&time[0],&res[0]);
  gdenoise->SetLineColor(kViolet-1);
  gdenoise->Draw("SAME L");


  TCanvas *c3 = new TCanvas();
  TGraph *gnoise_speech = new TGraph(n2,&time2[0],&speech[0]);
  TGraph *gfiltered_speech = new TGraph(n2,&time2[0],&res2[0]);
  TGraph *gcompare = new TGraph(n2,&speech[0],&res2[0]);
  gnoise_speech->SetLineColor(kRed);
  gfiltered_speech->SetLineColor(kBlue);
  gnoise_speech->Draw("AL");
  gfiltered_speech->Draw("SAME L");

  TCanvas *c4 = new TCanvas();
  gcompare->Draw("AP");
  
  ofstream fout;
  fout.open("filtered_speech.dat",ios::out);

  fout << "; Sample Rate 48000 \n";
  fout << "; Channels 1 \n";
  for(Int_t i =0; i<n2; i++){
    fout << time2[i] << " " << res2[i] << endl;
  }
  fout.close();


 
  
}
