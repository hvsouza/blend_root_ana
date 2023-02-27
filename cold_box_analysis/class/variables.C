// ________________________________________ //
// Author: Henrique Souza
// Filename: variables.C
// Created: 2021
// ________________________________________ //

template<Int_t npts = memorydepth>
class ADC_DATA{
  public:
    Double_t peak;
    Double_t peakpos;
    Double_t charge;
    Double_t fprompt;
    Double_t event;
    Double_t time;
    Double_t wvf[npts];
    Double_t base;
    Int_t selection;
    
    ADC_DATA();
    string tobranch = "peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[" + to_string(npts) + "]/D:base/D:selection/I";

    void setBranchName(Int_t n = npts){
      tobranch = "peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[" + to_string(n) + "]/D:base/D:selection/I";
    }

    // const char *tobranch = Form("peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[%i]/D:selection/I",npts);
};


template <Int_t npts>
ADC_DATA<npts>::ADC_DATA() {

  for(Int_t i = 0; i<npts; i++){
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
