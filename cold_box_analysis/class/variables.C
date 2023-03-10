// ________________________________________ //
// Author: Henrique Souza
// Filename: variables.C
// Created: 2021
// ________________________________________ //

class ADC_DATA {
  public:
    Double_t peak;
    Double_t peakpos;
    Double_t charge;
    Double_t fprompt;
    Double_t event;
    Double_t time;
    Double_t wvf[memorydepth];
    Double_t base;
    Int_t selection;

    ADC_DATA();
    string tobranch = "peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[" + to_string(memorydepth) + "]/D:base/D:selection/I";

    void setBranchName(Int_t n = memorydepth){
      tobranch = "peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[" + to_string(n) + "]/D:base/D:selection/I";
    }

    // const char *tobranch = Form("peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[%i]/D:selection/I",memorydepth);
};


ADC_DATA::ADC_DATA() {

  for(Int_t i = 0; i<memorydepth; i++){
    wvf[i]=0;
  }
  peak = 0;
  peakpos = 0;
  charge = 0;
  fprompt = 0;
  event = 0; 
  time = 0;
  selection = 0;
  base = 0;
}




class MY_DATA{
public:
  Int_t npts;
  Double_t peak;
  Double_t peakpos;
  Double_t charge;
  Double_t fprompt;
  Double_t event;
  Double_t time;
  Double_t *wvf;//[npts]
  Double_t base;
  Int_t selection;

  MY_DATA();
  ~MY_DATA();
  void Set_npts(Int_t n = 100);
};

MY_DATA::MY_DATA() {
  peak = 0;
  peakpos = 0;
  charge = 0;
  fprompt = 0;
  event = 0;
  time = 0;
  selection = 0;
  base = 0;
  npts = 0;
  wvf = 0;
}

MY_DATA::~MY_DATA() {
  delete [] wvf;
}

void MY_DATA::Set_npts(Int_t n) {
  if (n == npts) return; // no change needed
  if ((n < 1) || (n > npts)) { // if possible, "reuse" the current array
    delete [] wvf;
    wvf = ((n > 0) ? new Double_t[n] : 0);
  }
  npts = ((n > 0) ? n : 0);
}
