// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer.C
// Created: 2022
// ________________________________________ //
//
#include "MYCODES.h"



class ANALYZER{
  public:

    TFile *f = nullptr;
    TTree *t1 = nullptr;
    Int_t nentries = 0;
    vector<TBranch*> b;
    vector<ADC_DATA*> ch;
    Int_t nchannels = 1;
    vector<Int_t> channels = {1,2};
    vector<string> schannel;
    Double_t dtime = 4;
    string myname;
    string filename = "analyzed.root";
    TEventList *lev = nullptr;

    Int_t kch = 0;
    Int_t currentEvent = 0;
    DENOISE dn;
    WIENER *w = nullptr;
    string filter_type = "default";

    Int_t n_points;
    vector<vector<Double_t>> raw;
    vector<vector<Double_t>> wvf;
    vector<TH1D*> haverage;
    vector<TH1D*> hfft;
    TH1D *h = nullptr;
    Double_t *time = nullptr;
    Double_t *tempraw = nullptr;

    string plot_opt = "AL";
    TGraph *gwvf;
    string xlabel = "Time (ns)";
    string ylabel = "Amplitude (ADC Channels)";
    TH2D *hpers = nullptr;
    bool plotShowEvent = false;

    Double_t ymin = 0;
    Double_t ymax = 0;
    Double_t xmin = 0;
    Double_t xmax = 0;
    Double_t temp_charge = 0;
    Double_t temp_max = 0;
    Double_t temp_pos = 0;

    TCanvas *cpers = nullptr;


    // _________________________ Methods related to the class maintenance _________________________ //

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    void setAnalyzer(string m_filename = "analyzed.root"){
      filename = m_filename;
      if(f == nullptr) f = new TFile(filename.c_str(),"READ");
      if(t1 == nullptr) t1 = (TTree*)f->Get("t1");
      TList *lb = (TList*)t1->GetListOfBranches();
      if(lev==nullptr) lev = new TEventList(Form("lev_%s",myname.c_str()),Form("lev_%s",myname.c_str()));
      // branch title of new data format is equal to "ChX."
      // while for the old one, it is equal to `tobranch` string
      // So old branch will have brackets `[]`
      string branchtitle = lb->At(0)->GetTitle();
      int  bracketpos = branchtitle.find("[");
      Bool_t is_old_data = (bracketpos != -1) ? true : false;

      nchannels = lb->GetEntries();

      b.resize(nchannels);
      schannel.resize(nchannels);
      channels.resize(nchannels);
      ch.resize(nchannels);

      for (Int_t i = 0; i < nchannels; i++) {
        schannel[i] = lb->At(i)->GetName();
        channels[i] = schannel[i][2] - '0';
        ch[i] = new ADC_DATA();
        string temp_ch_name = schannel[i];
        b[i] = t1->GetBranch(temp_ch_name.c_str());
        b[i]->SetAddress(&ch[i]);
      }
      if (is_old_data){
        cout << "@@@@@@@@@@@@@@@ ________________________________________________________ @@@@@@@@@@@@@@@" << endl;
        cout << "@@@@@@@@@@@@@@@ This data format is not used anymore   " << endl;
        cout << "@@@@@@@@@@@@@@@ You can reprocess the data (no data should be loss) " << endl;
        cout << "@@@@@@@@@@@@@@@ A copy of the orignal file will be created as backup_" << filename << endl;
        cout << "@@@@@@@@@@@@@@@ Make sure that `memorydepth` is set correctly " << endl;
        cout << "@@@@@@@@@@@@@@@ ________________________________________________________ @@@@@@@@@@@@@@@" << endl;
        cout << "\n" << endl;

        string input;
        cout << "Do you want to recreate the root file? ([Y]es/[N]o) " << endl;
        cin >> input;
        if(input == "y" || input == "Y" || input == "Yes" || input == "yes"){
          // gets "peak/D:peakpos/D:charge/D:fprompt/D:event/D:time/D:wvf[2500]/D:base/D:selection/I"

          // Find the position of '[' and ']' in the string
          size_t left_bracket_pos = branchtitle.find("[");
          size_t right_bracket_pos = branchtitle.find("]");

          // Extract the substring between '[' and ']'
          string number_str = branchtitle.substr(left_bracket_pos + 1, right_bracket_pos - left_bracket_pos - 1);

          // Convert the extracted substring to an integer
          int number = std::stoi(number_str);
          if(number != memorydepth){
            cout << "\nERROR: the array length declared as `memorydepth` does not match the data array length" << endl;
            cout << "Declared: " << memorydepth << ", Actual size: " << number << endl;
            cout << "To solve this, quit root and execute again as:" << endl;
            cout << "myclass " << filename << " " << number << endl;
            return;
          }

          recreateFile();
        }
        else{
          cout << "Functions will not work" << endl;
          return;
        }

      }
      else{
        b[0]->GetEntry(0);
        n_points = ch[0]->npts;
      }
      setEmpty();
      nentries = t1->GetEntries();
    }

    void setEmpty(){
      string wienername = myname + "_w";
      if(!w) w = new WIENER(wienername.c_str(),dtime,250,1e-9,1e6,n_points);
      raw.resize(nchannels);
      wvf.resize(nchannels);
      haverage.resize(nchannels);
      hfft.resize(nchannels);
      for(Int_t i = 0; i < nchannels; i++){
        raw[i].resize(n_points);
        wvf[i].resize(n_points);
        hfft[i] = nullptr;
      }
      time = new Double_t[n_points];
      tempraw = new Double_t[n_points];
      for (int j = 0; j < n_points; j++) {
        time[j] = j*dtime;
      }
      xmin = 0;
      xmax = n_points*dtime;
    }

    Int_t getIdx(){
      return channels[kch];
    }
    void print(){
      t1->Print();
    }
    Bool_t setChannel(string mych = "Ch1."){
      for(Int_t i = 0; i < nchannels; i++){
        if(mych == schannel[i])
        {
          kch = i;
          return true;
        }
      }

      printf("%s not found, run s.print() to check the branches\n",mych.c_str());
      return false;
    }


    void setAnalyzerExt(TFile *&ft, TTree *&tr, vector<TBranch *> &bt){
      f = ft;
      tr = (TTree*)ft->Get("t1");
      t1 = tr;
      bt.resize(nchannels);
      setAnalyzer();
      for (Int_t i = 0; i < nchannels; i++) bt[i] = b[i];
    }

    void recreateFile();

    // _________________________ ________________________________________ _________________________ //


    // _________________________ Methods that get me something _________________________ //
    
    void getFFT(Double_t *_v = nullptr, bool inDecibel = false){
      if(_v == nullptr) _v = ch[kch]->wvf;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_v[i]);
      }
      w->fft(w->hwvf);
      if(inDecibel) w->convertDecibel();
      h = w->hfft;
    }

    Double_t getMean(Double_t from, Double_t to){
      if (from == to) return 0;
      Double_t res = 0;
      Double_t totalev = 0;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        res += ch[kch]->wvf[i];
        totalev+=1;
      }
      return res/totalev;
    }

    Double_t getMaximum(Double_t from, Double_t to, Double_t *v = nullptr){
      if (!v) v = ch[kch]->wvf;
      Double_t max = -1e12;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        if(v[i]>max){
          max = v[i];
          temp_max = max;
          temp_pos = i;
        }
      }
      return max;
    }

    void gen_rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, Double_t filter = 0, TH1D *htemp = nullptr, string selection = ""){
      getSelection(selection);
      kch = channel;
      htemp->GetYaxis()->SetTitle("# of events");
      htemp->GetXaxis()->SetTitle("Rise time_{90} (ns)");

      for(Int_t i = 0; i < lev->GetN(); i++){
        getWaveform(lev->GetEntry(i),kch);
        applyDenoise(filter);
        Double_t val = rise_time(kch, baseline_range, ispulse, peak_range);
        htemp->Fill(val);
      }
    }

    Double_t linear_interpole_tot(Int_t i, Double_t th){
      if (i - 1 < 0) {
        return i;
      }
      Double_t t1 = i-1;
      Double_t t2 = i;
      Double_t A1 = ch[kch]->wvf[i-1];
      Double_t A2 = ch[kch]->wvf[i];

      return ((t2-t1)*th - (A1*t2- A2*t1))/(A2-A1);
    }

    Double_t getTOT(Int_t channel = 0, Double_t from = 0, Double_t to = 0, Double_t time_trigger = 20, Double_t minimum_gap_ns = 12){
      kch = channel;
      Bool_t triggered = false;
      Double_t time_mark = 0;
      for(Int_t i = from/dtime; i<to/dtime; i++){
        if(ch[channel]->wvf[i]>time_trigger && triggered == false){
          triggered = true;
          time_mark = linear_interpole_tot(i, time_trigger)*dtime;
          temp_pos = time_mark;
          // break;
        }
        else if(triggered==true && ch[channel]->wvf[i]<time_trigger && (i*dtime-time_mark)>=minimum_gap_ns){
          time_mark = linear_interpole_tot(i, time_trigger)*dtime - time_mark;
          return time_mark;
        }
      }
      return 0;
    }

    void integrate(Double_t from = 0, Double_t to = 0, Double_t percent = 0, Bool_t withlimit = false){
      Double_t res = 0;
      Double_t max = -1e12;
      if (to == 0) to = n_points*dtime;
      if(percent == 0){ // normal integration
        for(Int_t i = from/dtime; i < to/dtime; i++){
          res += ch[kch]->wvf[i];
          if(ch[kch]->wvf[i]>=max){
            max = ch[kch]->wvf[i];
          }

        }
      }
      else{
        max = getMaximum(from, to);
        Int_t imin = 0;
        Int_t imax = n_points;
        if(withlimit){
          imin = from/dtime;
          imax = to/dtime;
        }
        for(Int_t i = temp_pos; i >= imin ;i--){
          Double_t val = ch[kch]->wvf[i];
          if(val >= percent*max){
            res += val;
          }
          else{
            break;
          }
        }
        for(Int_t i = temp_pos+1; i < imax ;i++){
          Double_t val = ch[kch]->wvf[i];
          if(val >= percent*max){
            res += val;
          }
          else{
            break;
          }
        }
      }
      temp_charge = res*dtime;
      temp_max = max;
    }
    void getIntegral(TH1D *_htemp = nullptr, Double_t from = 0, Double_t to = 0, string selection = "", Double_t filter = 0, Double_t sphe = 1., Double_t percent = 0){
      getSelection(selection);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        integrate(from, to, percent);
        _htemp->Fill(temp_charge/sphe);
      }
    }

    void exportAsHistogram(TH1D *hexport = nullptr, Double_t x_start = 0){
      TH1D *_htemp = new TH1D("","",n_points, x_start-dtime/2, x_start+n_points*dtime - dtime/2);
      for(Int_t i = 0; i < n_points; i++){
        _htemp->SetBinContent(i+1, ch[kch]->wvf[i]);
      }
      if(hexport==nullptr) h = _htemp;
      else hexport = _htemp;
    }

    void getWaveFromHistogram(TH1D *htemp){
      if (htemp->GetNbinsX() != n_points){
        cout << "Not same amount of samples! Graph has " << htemp->GetNbinsX() << endl;
        return;
      }
      for(Int_t i = 0; i < htemp->GetNbinsX(); i++){
        ch[kch]->wvf[i] = htemp->GetBinContent(i+1);
      }
    }
    void getWaveFromGraph(TGraph *gtemp){
      Double_t *xtemp = nullptr;
      Int_t ntemp = gtemp->GetN();
      if (ntemp != n_points){
        cout << "Not same amount of samples! Graph has " << ntemp << endl;
        return;
      }
      for(Int_t i = 0; i < ntemp; i++){
        ch[kch]->wvf[i] = *(gtemp->GetY()+i);
      }
    }


    void getWaveform(Int_t myevent = 0, Int_t k = -1){
      if (k == -1) k = kch;
      if (k>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return;
      }
      kch = k;
      b[k]->GetEvent(myevent);
      n_points = ch[k]->npts;
      currentEvent = ch[k]->event;
    }

    bool getWaveformHard(Int_t myevent = 0, Double_t factor = 1){
      if (kch>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return false;
      }
      b[kch]->GetEvent(myevent);
      n_points = ch[kch]->npts;
      currentEvent = ch[kch]->event;
      for (int j = 0; j < n_points; j++) {
        raw[kch][j] = ch[kch]->wvf[j]*factor;
        ch[kch]->wvf[j] = ch[kch]->wvf[j]*factor;
        wvf[kch][j] = raw[kch][j];
        time[j] = j*dtime;
      }
      return true;

    }

    void draw_rise_lines(Double_t start, Double_t finish, Double_t baseline_level, Double_t peak_level){
      Double_t gxmin = gPad->GetUxmin();
      Double_t gxmax = gPad->GetUxmax();
      TLine *lbase = new TLine(gxmin, baseline_level, gxmax, baseline_level);
      TLine *lpeak = new TLine(gxmin, peak_level, gxmax, peak_level);
      TLine *l90 = new TLine(gxmin, 0.9*(peak_level-baseline_level)+baseline_level, gxmax, 0.9*(peak_level-baseline_level)+baseline_level);
      TLine *l10 = new TLine(gxmin, 0.1*(peak_level-baseline_level)+baseline_level, gxmax, 0.1*(peak_level-baseline_level)+baseline_level);

      TLine *lstart = new TLine(start, baseline_level, start, peak_level);
      TLine *lfinish = new TLine(finish, baseline_level, finish, peak_level);

      set_up_lines(lbase,kGreen+2);
      set_up_lines(lpeak,kGreen+2);
      set_up_lines(l90, kRed);
      set_up_lines(l10, kRed);
      set_up_lines(lstart, kRed);
      set_up_lines(lfinish, kRed);

    }

    Double_t triggerTime(vector<Double_t> baseline_range = {0,0}, vector<Double_t> peak_range = {0,0}, Double_t threshold = 0, bool debug = false){
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);

      Int_t startpoint = peak_range[0]/dtime;
      Int_t endpoint = peak_range[1]/dtime;

      Double_t time_mark;
      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch]->wvf[i] >= threshold) {
          if(i == startpoint) break;
          time_mark = linear_interpole_tot(i, threshold - baseline_level)*dtime;
          temp_pos = time_mark;
          return time_mark;
        }
      }
      temp_pos = 0;
      return 0;
    }


    Double_t rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, bool debug = false){
      kch = channel;
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);

      Double_t peak_level = 0;
      if (ispulse){
        peak_level = getMaximum(peak_range[0],peak_range[1]);
      }
      else{
        peak_level = getMean(peak_range[0],peak_range[1]);
      }
      Int_t startpoint = peak_range[0]/dtime;
      Int_t endpoint = peak_range[1]/dtime;
      Bool_t triggered = false;
      Double_t time_mark = 0;

      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch]->wvf[i] >= 0.1*(peak_level-baseline_level)+baseline_level && triggered==false) {
          triggered = true;
          time_mark = linear_interpole_tot(i, 0.1*(peak_level-baseline_level))*dtime;
          // time_mark = i*dtime;
        }
        else if(triggered == true && ch[kch]->wvf[i]>=0.9*(peak_level - baseline_level)+baseline_level){
          if(debug) draw_rise_lines(time_mark, i*dtime, baseline_level, peak_level);
          temp_pos = linear_interpole_tot(i, 0.9*(peak_level - baseline_level))*dtime;
          // time_mark = i*dtime - time_mark;
          temp_max = peak_level;
          return temp_pos - time_mark;
        }
      }
      if(debug) draw_rise_lines(time_mark, 0, baseline_level, peak_level);
      return 0;
    }

    void getBackFFT(Double_t *_filtered = nullptr){
      if (_filtered == nullptr){_filtered = ch[kch]->wvf;}
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }

    }
    void averageFFT(Int_t maxevent = 0, string selection = "", bool inDecibel = true, Double_t factor = 0, Double_t filter = 0){
      if (maxevent==0) {
        maxevent = nentries;
      }
      getSelection(selection);
      Int_t nev = lev->GetN();
      if (maxevent < nev) {
        nev = maxevent;
      }
      Int_t iev = 0;
      hfft[kch] = (TH1D*)w->hfft->Clone(Form("hfft_%s_ch%d",myname.c_str(),kch));
      hfft[kch]->Reset();
      cout << "\n";
      Int_t total = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        if(i%200==0) cout << "computing event " << i << " of " << nev << "\r" << flush;
        getWaveform(iev,kch);
        applyDenoise(filter);
        // applyFreqFilter();
        getFFT();
        for (Int_t j = 0; j < n_points/2; j++) hfft[kch]->AddBinContent(j+1,w->hfft->GetBinContent(j+1));
        total++;
      }
      cout << "\n";
      hfft[kch]->Scale(1./total);
      hfft[kch]->SetEntries(total);

      if (inDecibel){
        w->convertDecibel(hfft[kch], factor);
        hfft[kch]->GetYaxis()->SetTitle("Magnitude (dB)");
      }
      h = hfft[kch];
    }

    void averageWaveform(Int_t maxevent = 0, string selection = "", Double_t filter =  0){
      if (maxevent==0) {
        maxevent = nentries;
      }
      if(!haverage[kch]) haverage[kch] = new TH1D(Form("haverage_%s_Ch%d",myname.c_str(),channels[kch]),"Averaged waveform",n_points,0,n_points*dtime);
      else haverage[kch]->Reset();
      haverage[kch]->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
      haverage[kch]->GetXaxis()->SetTitle("Time (ns)");
      Int_t total = 0;
      getSelection(selection);
      Int_t nev = lev->GetN();
      if (maxevent < nev) {
        nev = maxevent;
      }
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        total += 1;
        for (Int_t j = 0; j < n_points; j++){
          haverage[kch]->AddBinContent(j+1,ch[kch]->wvf[j]);
        }

      }
      haverage[kch]->Scale(1./total);
      haverage[kch]->SetEntries(total);
      haverage[kch]->Sumw2(kFALSE);
      h = haverage[kch];
    }

    void zeroCrossSearch(Double_t *derwvf, vector<Int_t> &peaksRise, vector<Int_t> &peaksCross, Int_t start, Int_t finish){
      if(start == 0) start = dtime;
      bool lastIsPositive = false; // I want to always start with a positive crossing
      for(Int_t i = start/dtime; i < finish/dtime-1; i++){
        if(lastIsPositive==false && derwvf[i] >= 0 && derwvf[i-1] <=0 && derwvf[i+1]>0){
          peaksRise.push_back(i);
          lastIsPositive = true;
          i = i+1;
        }
        else if(lastIsPositive && derwvf[i] <= 0 && derwvf[i-1] >=0 && derwvf[i+1] < 0){
          peaksCross.push_back(i);
          lastIsPositive = false;
        }
      }
    }


    void convert_1D_histogram(TH2D *horiginal, TH1D **hnew, vector<vector<Double_t>> &_hy, vector<vector<Double_t>> &_hx, Int_t rebin = 1){
      Int_t nx = horiginal->GetNbinsX();
      Int_t ny = horiginal->GetNbinsY();
      Double_t hmin = horiginal->GetXaxis()->GetBinLowEdge(1);
      Double_t hmax = horiginal->GetXaxis()->GetBinLowEdge(nx) + horiginal->GetXaxis()->GetBinWidth(nx);

      *hnew = new TH1D("hnew","hnew",nx/rebin,hmin,hmax);
      Double_t sum = 0;
      Double_t div = 0;
      Double_t stddev = 0;
      Double_t n = 0;
      Int_t nrebins = 0;
      int ndiv = 0;
      for(Int_t i = 0; i < nx; i++) { // for loop in the X axis of the histogram;
        for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
          Double_t weight = horiginal->GetBinContent(i+1,j+1);
          if (weight != 0){
            Double_t bnpes = horiginal->GetYaxis()->GetBinCenter(j);
            sum += weight*bnpes;
            div += weight;
          }
        }
        ndiv++;
        if(ndiv == rebin){
          if(div == 0){sum = 0; div = 1;}
          for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
            Double_t weight = horiginal->GetBinContent(i+1,j+1);
            if (weight != 0){
              Double_t bnpes = horiginal->GetYaxis()->GetBinCenter(j);
              stddev += weight*(bnpes-sum/div)*(bnpes-sum/div);
              n += 1;
            }
          }
          if (n == 1){ n = 2;}
          stddev = (stddev/div)*n/(n-1);
          (*hnew)->SetBinContent(nrebins+1, sum/div);

          Double_t err = sqrt(abs(sum)/div);
          (*hnew)->SetBinError(nrebins+1,sqrt(stddev + err*err));
          if(sum!=0){
            _hy[0].push_back(sum/div);
            _hy[1].push_back((*hnew)->GetBinError(nrebins+1));
            _hx[0].push_back((*hnew)->GetBinCenter(nrebins+1));
            _hx[1].push_back((*hnew)->GetBinWidth(nrebins+1)/2.);
          }

          nrebins++;

          sum = 0;
          div = 0;
          stddev = 0;
          n = 0;
          ndiv = 0;
        }


      }

    }
    // _________________________ _____________________________ _________________________ //


    // _________________________ Filters and processing _________________________ //
    void scaleWvf(Double_t factor, Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      for (Int_t i = 0; i < n_points; i++) {
        _filtered[i] = _filtered[i]*factor;
      }
    }

    void addOffet(Double_t offset = 0, Double_t from = 0, Double_t to = 0){
      if(offset == 0){
        offset = ch[kch]->base;
      }
      if(to == 0) to = n_points*dtime;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        ch[kch]->wvf[i] = ch[kch]->wvf[i] + offset;
      }
    }
    void checkSignals(Double_t **_raw, Double_t **_filtered){
      if (*_raw == nullptr  && *_filtered == nullptr){
        *_filtered = ch[kch]->wvf;
        *_raw = tempraw;
          for (Int_t i = 0; i < n_points; i++) {
            (*_raw)[i] = ch[kch]->wvf[i];
          }
      }
      else if(*_raw == *_filtered){
        for (Int_t i = 0; i < n_points; i++) {
          tempraw[i] = (*_raw)[i];
        }
        *_raw = tempraw;
      }
    }
    void differenciate(Double_t factor = 1, Double_t *_raw = nullptr, Double_t *_shifted = nullptr){
      checkSignals(&_raw,&_shifted);
      _shifted[0] = 0;
      _shifted[n_points-1] = 0;
      for(int i=1; i<n_points-1; i++){
        _shifted[i] = factor*(_raw[i+1] - _raw[i-1])/(2*dtime);
        // _shifted[i]=_raw[i] - (i-delay_time>=0 ? _raw[i-delay_time] : 0);
      }
    }

    void applyMovingAverage(Int_t mafilter = 0, Double_t start = 0, Double_t finish = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      if (finish == 0) finish = n_points*dtime;
      if(mafilter!=0) {
        checkSignals(&_raw,&_filtered);
        dn.movingAverage(_raw,_filtered,mafilter,start/dtime,finish/dtime);
      }
    }

    void applyDenoise(Double_t filter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      checkSignals(&_raw,&_filtered);
      if (filter == 0 && _filtered == ch[kch]->wvf) return;
      if (filter_type == "default"){
        applyTV1D(filter, _raw, _filtered);
      }
      else if (filter_type == "ma"){
        applyMovingAverage(filter);
      }
      else if(filter_type == "band")
        applyBandCut();
      else{
        applyFreqFilter();
      }
    }

    void applyTV1D(Double_t filter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      checkSignals(&_raw,&_filtered);
      if (filter == 0) return;
      dn.TV1D_denoise(_raw,_filtered,n_points,filter);
    }

    void setFreqFilter(Double_t frequency_cut, string filter_type = "gaus"){
      this->filter_type = filter_type;
      w->setFilter(frequency_cut,filter_type);
    }
    void applyFreqFilter(Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->apply_filter();
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
      // w->hwvf->Draw("");
    }


    void applyBandCut(WIENER *_temp = nullptr, Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      if(_temp == nullptr) _temp = w;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->applyBandCut(_temp);
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
    }

    void makeCopy(Double_t *cpy, Double_t *original = nullptr){
      if (original == nullptr) original = ch[kch]->wvf;
      for (Int_t i = 0; i < n_points; i++) {
        cpy[i] = original[i];
      }
    }

    void shift_waveform(TH1D *h, Int_t new_max, Bool_t rawShift = false){
      Int_t old_max = h->GetMaximumBin();
      if(rawShift) old_max = 0;
      Int_t old_ref = old_max - new_max;
      TH1D *htemp = (TH1D*)h->Clone("htemp");
      Double_t temp;
      if(old_ref<0){
        // cout << " case lower" << endl;
        old_ref = n_points-(new_max-old_max);
      }
      for(Int_t i = 1; i<n_points-(old_ref); i++){
        temp = htemp->GetBinContent(old_ref+i);
        h->SetBinContent(i,temp);
      }
      Int_t aux = 1;
      for(Int_t i = n_points-(old_ref); i<=n_points; i++){
        temp = htemp->GetBinContent(aux);
        h->SetBinContent(i,temp);
        aux++;
      }
      delete htemp;
    }


    Double_t reval_baseline(vector<Double_t> range_base, Double_t exclusion_baseline, Double_t exclusion_window, TH1D *hbase = nullptr){
      Double_t result = 0;
      if(!hbase) hbase = new TH1D("hbase","finding baseline",TMath::Power(2,14),0,TMath::Power(2,14));
      hbase->Reset();
      for(Int_t i=range_base[0]/dtime; i<range_base[1]/dtime; i++) hbase->Fill(ch[kch]->wvf[i]);
      Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());
      Double_t hmean = hbase->GetMean();
      Double_t hstd = hbase->GetStdDev();
      Double_t bins=0;
      for(Int_t i=range_base[0]/dtime; i<range_base[1]/dtime;){
        if(ch[kch]->wvf[i] > res0 + exclusion_baseline || ch[kch]->wvf[i]<res0 - exclusion_baseline) {
          i+=exclusion_window/dtime;
        }
        else{
          result += ch[kch]->wvf[i];
          bins+=1;
          i++;
        }
      }
      if(bins>0)result/=bins;
      if(bins > (range_base[1]/dtime)/3.){
      }
      else{
        result = res0;
        cout << "Not enough points probably.. " << endl;
      }
      addOffet(-result);
      return result;
    }



    // _________________________ ______________________ _________________________ //

    // _________________________ Methods for selection _________________________ //
    bool checkHigher(Double_t a, Double_t b){
      if (a > b){return true;}
      else{return false;}
    }
    bool checkLower(Double_t a, Double_t b){
      if (a < b){return true;}
      else{return false;}
    }

    void invertSelection(){
      TEventList *ttemp = new TEventList("ttemp", "ttemp");
      t1->Draw(">>ttemp","");
      ttemp->Subtract(lev);
      lev->Reset();
      for (Int_t i = 0; i < ttemp->GetN(); i++) {
        lev->Enter(ttemp->GetEntry(i));
      }
      delete ttemp;
    }

    void selectByAmplitude(Double_t filter = 0, Double_t xmin = 0, Double_t xmax = 0, Double_t limit = 100, string type = "higher"){
      f->cd();
      lev = (TEventList*)gDirectory->Get(Form("lev_%s",myname.c_str()));
      if(lev->GetN() == 0){
        getSelection("");
      }
      TEventList *ttemp = new TEventList("ttemp", "ttemp");

      if (xmax == 0) {
        xmax = n_points*dtime;
      }
      for(Int_t i = 0; i < lev->GetN(); i++){
        Int_t ev = lev->GetEntry(i);
        getWaveform(ev, kch);
        applyDenoise(filter);
        bool contact = false;
        for(Int_t j = xmin/dtime; j < xmax/dtime; j++){
          if(type == "higher") contact = checkHigher(ch[kch]->wvf[j], limit);
          else contact = checkLower(ch[kch]->wvf[j], limit);
          if(contact){
            ttemp->Enter(ev);
            break;
          }
        }
      }
      lev->Subtract(ttemp);
      delete ttemp;
    }

    void getSelection(string selection)
    {
      f->cd();
      if(selection!="use_mine"){
        t1->Draw(Form(">>lev_%s",myname.c_str()),selection.c_str());
      }
      else{
        if (lev->GetN() == 0){
          cout << "No event selectected... " << endl;
        }
      }
    }
    // _________________________ _____________________ _________________________ //

    void set_up_lines(TLine *l, Color_t color){
      l->SetLineColorAlpha(color,0.6);
      l->SetLineStyle(9);
      l->Draw("SAME");

    }


    Double_t computeSNR_simple(Double_t xmin, Double_t xmax, vector<Double_t> signal_range){
      Double_t snr = 0;
      Double_t avg = 0;
      Double_t sum = 0;
      Double_t navg = 0;
      Double_t stddev = 0;

      for(Int_t i = xmin/dtime; i < xmax/dtime; i++){
        if(i*dtime > signal_range[1] || i*dtime < signal_range[0]){
          sum += ch[kch]->wvf[i];
          navg += 1;
        }
      }
      avg = sum/navg;
      for(Int_t i = xmin/dtime; i < xmax/dtime; i++){
        if(i*dtime > signal_range[1] || i*dtime < signal_range[0]){
          Double_t diff = (ch[kch]->wvf[i] - avg);
          stddev += diff*diff;
        }
      }
      stddev = sqrt(stddev/(navg-1));
      Double_t max = getMaximum(signal_range[0], signal_range[1]);
      snr = max/stddev;
      return snr;
    }




    void showFFT(Int_t naverage = 10, Int_t maxevent = 0, Int_t dt = 100, bool inDecibel = true);
    void sample_plot(Int_t myevent = 0, string opt = "", Double_t filter = 0, Double_t factor = 1., Int_t mafilter = 0);
    void showWaveform(Int_t maxevent = 0, Double_t filter = 0, Int_t dt = 0);
    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Double_t filter = 0, string cut="", Double_t factor = 1);
    void add_persistence_plot(TH2D *_htemp = nullptr, Double_t filter = 0, string cut = "", Double_t factor = 1);
    TGraph drawGraph(string opt = "", Int_t n = 0, Double_t* x = nullptr, Double_t* y = nullptr);
    void debugSPE(Int_t event, Int_t moving_average, Int_t n_moving, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> not_used, Double_t filter = 16, Double_t factor = 200, Double_t *SNRs = nullptr);
    void minimizeParamsSPE(Int_t event, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> rangeInter = {0,0}, Double_t filter = 16, Double_t factor = 200);
    void drawZeroCrossingLines(vector<Int_t> &peaksCross, vector<Int_t> &peaksRise, TCanvas *c = nullptr, Double_t ymin = 0, Double_t ymax = 0);
    void histoTimeTrigger(Int_t nstart = 0, Int_t nfinish = 0, TH1D *_htemp = nullptr);
    void graphTimeTrigger(Int_t nstart = 0, Int_t nfinish = 0, TGraph *_gtemp = nullptr);
    void check_filtering(vector<Int_t> filter_max_and_step = {0,0}, Int_t event = 0, Int_t rebine = 1);

    ANALYZER(string m_myname = "z") : myname{m_myname}{
    }

};



void ANALYZER::recreateFile(){
  TFile *ftemp = new TFile("new_data.root","RECREATE");
  TTree *ttemp = new TTree("t1","ADC processed waveform");

  vector<TBranch*> bold(nchannels);
  vector<ADC_DATA*> newch(nchannels);
  vector<OLD_ADC_DATA*> oldch(nchannels);



  n_points = memorydepth;
  TFile *oldf = new TFile(filename.c_str(),"READ");
  TTree *oldt = (TTree*)oldf->Get("t1");

  for(Int_t i = 0; i < (int)channels.size(); i++){
    oldch[i] = new OLD_ADC_DATA();
    ttemp->Branch(Form("Ch%i.",channels[i]),&newch[i]);
    bold[i] = oldt->GetBranch(Form("Ch%i",channels[i]));
    bold[i]->SetAddress(oldch[i]);

  }

  cout << "\n";
  Int_t nevts = oldt->GetEntries();
  for(Int_t j = 0; j < nevts; j++){
    if(j%200==0) cout << "computing event " << j << " of " << nevts << "\r" << flush;
    for(Int_t i = 0; i < (int)channels.size(); i++){
      bold[i]->GetEntry(j);

      newch[i]->peak      = oldch[i]->peak;
      newch[i]->peakpos   = oldch[i]->peakpos;
      newch[i]->charge    = oldch[i]->charge;
      newch[i]->fprompt   = oldch[i]->fprompt;
      newch[i]->event     = oldch[i]->event;
      newch[i]->time      = oldch[i]->time;
      newch[i]->selection = oldch[i]->selection;
      newch[i]->base      = oldch[i]->base;
      newch[i]->Set_npts(n_points);
      for(Int_t k = 0; k < n_points; k++){
        newch[i]->wvf[k] = oldch[i]->wvf[k];

      }
    }

    ttemp->Fill();
  }
  cout << "\n";

  ftemp->WriteTObject(ttemp,"t1","TObject::kOverwrite");

  TIter next(oldf->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    string lname = key->GetName();
    if(lname == "t1") continue;
    auto *tobj = oldf->Get(lname.c_str());
    ftemp->WriteTObject(tobj, lname.c_str());
  }

  delete ttemp;
  delete oldt;
  ftemp->Close();
  oldf->Close();
  delete ftemp;
  delete oldf;
  for(Int_t i = 0; i < (int)channels.size(); i++){
    delete newch[i];
    delete oldch[i];
  }
  f->Close();
  string make_backup = Form("mv %s backup_%s",filename.c_str(),filename.c_str());
  string rename = Form("mv new_data.root %s",filename.c_str());

  int outputsys = system(make_backup.c_str());
  if(outputsys != 0){
    cout << "\n\nCould not backup properly. New data is saved as \"new_data.root\"" << endl;
    filename = "new_data.root";
  }
  else{
    outputsys = system(rename.c_str());
  }

  f = nullptr;
  t1 = nullptr;
  lev = nullptr;
  setAnalyzer(filename.c_str());


}
