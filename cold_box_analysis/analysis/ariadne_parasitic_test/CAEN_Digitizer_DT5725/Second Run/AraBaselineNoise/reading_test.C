#include <fstream>
#include <iostream>

void reading_test(){
  int headbin;
  int nbytes = 4; 
  int memorydepth = 0;
  uint32_t valbin = 0;
  vector<Double_t> raw;
  Bool_t first_line = true;
  
  ifstream fin;
  fin.open("0_wave0.dat", ios::in | ios::binary);

  if(fin.good() && fin.is_open()){ // Ok
    cout << "Reading file" << endl;
  }
  else{ 
    cout << "File did not open!!" << endl;
    return;      
  }
  while(!fin.fail()){
    for(Int_t ln=0;ln<6;ln++){ // 4 bytes (32 bits) for each head 
      fin.read((char *) &headbin, nbytes);
      // header0 will be EventSize, so: you can do
      if(ln==0){
        memorydepth = headbin;
        // the result is in bytes, so 2*NSamples+4*headers
        memorydepth = (memorydepth-4*6)/2;
        printf("File size: %d \n",memorydepth);
      }
    }
    if(first_line){
      raw.resize(memorydepth);
      first_line=false;
    }
    for(int j = 0; j < memorydepth; j++)
      {
        fin.read((char *) &valbin, 2); // 2 bytes (16 bits) per sample
        if(fin.bad() || fin.fail()){
          break;
        }
        raw[j] = valbin;
        printf("%d %.0f \n",j,raw[j]);
      }
  }

  
}
