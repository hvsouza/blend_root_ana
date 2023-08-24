// ________________________________________ //
// Author: Henrique Souza
// Filename: variables.C
// Created: 2021
// ________________________________________ //

#ifndef ADC_DATA_H
#define ADC_DATA_H

class ADC_DATA{
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

    ADC_DATA();
    virtual ~ADC_DATA();
    void Set_npts(Int_t n = 100);

};

ADC_DATA::ADC_DATA() {
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

ADC_DATA::~ADC_DATA() {
  delete [] wvf;
}

void ADC_DATA::Set_npts(Int_t n) {
  if (n == npts) return; // no change needed
  if ((n < 1) || (n > npts)) { // if possible, "reuse" the current array
    delete [] wvf;
    wvf = ((n > 0) ? new Double_t[n] : 0);
  }
  npts = ((n > 0) ? n : 0);
}


#endif
