void convert_data(){


  ifstream fin;
  ofstream fout; 
  fin.open("july8_led6p8_ch2_36V_1.txt",ios::in);
  // fin.open("july8_led6p8_ch3_36V_1.txt",ios::in);
  fout.open("working/july8_ch2.txt",ios::out);
  // fout.open("working/july8_ch3.txt",ios::out);

  if(!(fin.good()  && fin.is_open())) cout << "filed to open file" << endl;


  Int_t nsize = 2000;
  vector<Double_t> wled(nsize);
  vector<Double_t> wvf(nsize);
  vector<Double_t> time(nsize);
  vector<Double_t> temporary;
  
  for(Int_t i = 0; i<nsize; i++){
    wvf[i] = 0;
    time[i] = i*5;
  }
  Double_t led, signal;
  Int_t aux_pts = 0;
  while (!fin.eof()){
    fin >> signal >> led;
    aux_pts++;
    temporary.push_back(signal);
    if(aux_pts<1000) continue;
    if(led > 1.0){
      Int_t currentsize = temporary.size();
      for(Int_t i = 0; i<nsize/2; i++){
        wvf[nsize/2-1-i] = temporary[currentsize - 1 - i];
      }
      for(Int_t i = 0; i<nsize/2; i++){
        fout << wvf[i] << "\n";
      }
      for(Int_t i = nsize/2; i<nsize; i++){
        fin >> signal >> led;
        wvf[i] = signal;
        fout << wvf[i] << "\n";
      }
      temporary.clear();
    }
  }
  
}
