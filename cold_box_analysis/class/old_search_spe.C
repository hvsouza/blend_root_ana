
class SPHE{




  public:

    TFile *fout=nullptr;

    TTree *tout = nullptr;

    Bool_t creation = true; //to verify creation of ttree branches


    TFile *fwvf = nullptr;
    TTree *twvf = nullptr;


    Double_t value = 0;
    Double_t desv = 0;
    Int_t channel = 1;
    ADC_DATA<memorydepth> ch;
    ADC_DATA<memorydepth> sample;


    // ____________________ Variables to calculate ____________________ //

    TH1D *hbase = new TH1D("hbase","histogram for baseline",5*800,-400,400);
    TH1D *hbase_smooth = new TH1D("hbase_smooth","histogram for baseline smoothed",5*800,-400,400);
    // TH1D *hcharge = new TH1D("hcharge","",100000,-50000,50000);
    TH1D *hcharge = new TH1D("hcharge","",50000,0,0);
    TH1D *hzero = new TH1D("hzero","",120000,-200000,2*1300000);
    TH1D *hnobase = new TH1D("hnobase","",120000,-200000,2*1300000);

    TH1D *hstat = new TH1D("hstat","",120000,-200000,2*1300000);


    Double_t dtime = 2.;

    Double_t mean = 0;
    Double_t stddev = 0;
    Double_t tolerance; // n sigmas
    Double_t threshold;
    Bool_t fix_threshold = false;
    Double_t baseLimit;

    Double_t filter = 0;
    Bool_t withfilter = true;

    Double_t timeLimit; // time after LED signal
    Double_t timeLow ; // integration time before peak
    Double_t timeHigh ; // integration time after peak

    Double_t start = 0;
    Double_t finish = memorydepth*dtime;

    Double_t noiseLow;
    Double_t noiseHigh;

    Double_t forBaseline = 2000;
    Double_t baseCharge = 0;

    Double_t social_distance = 2;
    Bool_t check_selection = true;

    Int_t interactions; // for moving avarage
    Int_t midpoint = 0; // midpoint that takes the avarage
    Int_t width = 0;
    Double_t sum=0;

    vector<Double_t> peakPosition;
    vector<Double_t> peakMax;

    Double_t charge = 0;
    Double_t treeCharge = 0;
    Double_t treePeak = 0;
    Double_t strikes = 0;
    Double_t ptsLow = 0;
    Double_t ptsHigh = 0;
    Bool_t discard = false;

    TH1D *hfilter = new TH1D("hfilter","",5*800,-400,400);
    Double_t devFilter;


    Int_t cut = 0;

    Int_t gap = 2000; // at the and, we wont look at the final 4000 ns
    Double_t lowerThreshold = -9999;
    Double_t maxHits = 15;
    Double_t nSigmas_integration = 0.5;
    Int_t atleast = 120; //at least this time in the first integration and second integration
    Int_t lowPtsLimit = 20;
    Bool_t just_a_test = false; //if true, will only look 5000 events;
    Bool_t led_calibration = false; // for led integration only

    Int_t just_this = 5000;

    // ____________________ matching filter ____________________

    vector<Double_t> mySample;
    vector<Double_t> matched;
    vector<Double_t> matchingAreas;

    Double_t shift = 0;
    Double_t matched_value = 0;
    Double_t higherValue = 0;
    Double_t tempPosition;


    Double_t area_off = 0;
    Double_t matchTrigger = 1.;

    Double_t from = 0;
    Double_t to = 0;


    Bool_t matching = false;
    Bool_t badBaseline = false; // this need to be used in the first Gamma run...

    // ____________________ Variables to show events ____________________ //

    const static Int_t nshow = 100;
    Int_t my_events[nshow];
    Double_t eventNow = 0;
    Int_t aux_events = 0;


    TGraph *g_smooth = nullptr;
    TGraph *g_normal = nullptr;
    TGraph *g_points = nullptr;

    vector<Double_t> temp_peak;
    vector<Double_t> peak_smooth;
    vector<Double_t> timeg;

    vector<Double_t> selected_peaks;
    vector<Double_t> selected_time;
    Int_t teste = 0;


    Bool_t fillg = true;

    Int_t mark = 0;

    Double_t xi = 0; // time interval for sampling (ns)
    Double_t xf = 1000;
    Double_t center = 300;
    Int_t nsample = (xf-xi)/2+1;
    Double_t counter = 0;

    TH2D *hsample = new TH2D("sample","peak vs time", nsample, xi, xf, 3000,-1500,1500);

    vector<Double_t> ymean;
    vector<Double_t> ynormal;
    vector<Double_t> xmean;

    // variables to find mean waveform
    Bool_t get_wave_form = false;
    vector<Double_t> mean_waveforms;
    Double_t wvfcharge;
    Bool_t valid;
    Int_t naverages;
    Double_t mean_before = 100;
    Double_t mean_after = 400;
    Int_t low_cut = 800; // to select sphe waveforms, inside this region, peak should not be bigger then val_cut
    //it was 920 before
    Int_t high_cut = 2000;
    Double_t val_cut = 6;
    TMultiGraph *gm = new TMultiGraph();
    TGraph *gwaveforms = nullptr;
    Double_t shifter = 12;

    Double_t sphe_charge;
    Double_t sphe_charge2;
    Double_t delta;
    Double_t deltaplus = 1.4;
    Double_t deltaminus = 1.2;
    Double_t sphe_std;

    Double_t sphe_charge_ch0;

    Double_t sphe_charge2_ch0;

    Double_t sphe_std_ch0;
    Double_t sphe_std_ch1;

    string charge_status = "";

    Bool_t darkNoise = false;

    Double_t baselineTime=2000;
    Double_t too_big = 1000;
    Double_t waiting = 0;
    Double_t pre_filter = 0;
    Int_t peakPerEvent = 0;
    // ____________________________________________________________________________________________________ //

    void giveMeSphe_darkCount(string name){
      temp_peak.resize(memorydepth);
      gROOT->SetBatch(kTRUE);
      fout = new TFile(Form("sphe_histograms_darkCount_Ch%i.root",channel),"RECREATE");
      tout = new TTree("t1","baseline info");
      darkNoise = true;
      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
      }
      twvf = new TTree("t1","mean waveforms");

      if(led_calibration==false){
        makeHistogram(name);
        fout->WriteObject(tout,"t1","TObject::kOverwrite");
        // fout->Close();
      }
      else{
        makeSimpleHistogram(name);
        // get_wave_form = false;
      }
      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",channel));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
      }
      channel=channel + 1;
    }


    // ____________________________________________________________________________________________________ //
    void giveMeSphe(string name){
      gROOT->SetBatch(kTRUE);

      fout = new TFile("sphe_histograms.root","RECREATE");
      tout = new TTree("t1","baseline info");
      darkNoise = false;
      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
      }
      twvf = new TTree("t1","mean waveforms");

      makeHistogram(name);

      fout->WriteObject(tout,"t1","TObject::kOverwrite");
      fout->Close();
      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",channel));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
      }
      channel=channel + 1;
    }

    // ____________________________________________________________________________________________________ //
    void integrateSignal(){
      charge = 0;
      Int_t npeaks = (int)peakPosition.size();
      value = 0;
      desv = 0;

      Int_t discardedPeaks = 0;

      Int_t aux_sample = 0;


      sphe_charge = sphe_charge_ch0; // wave0
      sphe_charge2 = sphe_charge2_ch0; // wave0
      delta = sphe_charge2 - sphe_charge;
      sphe_std = sphe_std_ch0;

      Double_t statcharge[npeaks];
      Double_t statpeak[npeaks];
      Double_t trackLow[npeaks];
      Double_t trackHigh[npeaks];
      Double_t trackStrikes[npeaks];
      Double_t statpos[npeaks];
      Bool_t discard_this[npeaks];

      vector<vector<Double_t>> temp_waveforms(npeaks);
      vector<Double_t> waveforms(memorydepth);
      //     vector<TGraph*> gwaveforms(npeaks);
      vector<Bool_t> notAGoodWaveform(npeaks);


      Int_t auxstat = 0;
      Bool_t notGood = false;
      Double_t lastTime = 0 ;
      Double_t peakStd = 0;

      Double_t thismanypointsLow = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t thismanypointsHigh = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t thismanypointsBase = 0; //this is to check the amount of points, noise may cause low points that is terrible
      Double_t totalmany = 0; //this helps speacially for matching...
      Bool_t next_is_bad = false;
      Double_t waitingInterval = -1;
      for(Int_t i = 0; i<npeaks; i++){
        statcharge[i] = 0;
        statpeak[i] = 0;
        trackStrikes[i] = 0;
        discard_this[i] = false;
        notGood = false;
        thismanypointsLow = 0;
        thismanypointsHigh = 0;
        thismanypointsBase = 0;
        totalmany = 0;


        if((peakPosition[i]+social_distance*timeHigh>=peakPosition[i+1] && i+1!=npeaks) || next_is_bad || peakPosition[i]<waitingInterval){

          selected_peaks.push_back(temp_peak.at(peakPosition[i]/dtime));
          selected_time.push_back(peakPosition[i]);
          discard_this[i] = true;
          for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
            if(temp_peak.at(j)>=statpeak[i]){
              statpeak[i] = temp_peak.at(j);
            }
          }
          if(statpeak[i]>too_big){

            waitingInterval = peakPosition[i]+waiting;
          }
          if(peakPosition[i]<waitingInterval && waitingInterval!=-1) next_is_bad = true;
          else next_is_bad = false;
          continue;
        }
        else{
          waitingInterval = -1;
          next_is_bad = false;
        }



        notAGoodWaveform[i] = false;



        Int_t strikes = 0;
        // integration for the back part
        for(Int_t j = (peakPosition.at(i))/dtime; j>= (peakPosition[i]-timeLow)/dtime; j--){
          charge += temp_peak.at(j);
          if(withfilter) statcharge[i]+= temp_peak.at(j);
          else statcharge[i]+= ch.wvf[j];
          if(temp_peak.at(j)>=statpeak[i]){
            statpeak[i] = temp_peak.at(j);
          }
          selected_peaks.push_back(temp_peak.at(j));
          selected_time.push_back(timeg.at(j));
          thismanypointsLow++;
          totalmany++;
          if(temp_peak[j]<lowerThreshold){
            strikes++;
          }
        }


        // integration for the front part
        for(Int_t j = (peakPosition.at(i))/dtime + 1; j<=(peakPosition[i]+timeHigh)/dtime; j++){
          charge += temp_peak.at(j);

          if(withfilter) statcharge[i]+= temp_peak.at(j);
          else statcharge[i]+= ch.wvf[j];
          if(temp_peak.at(j)>=statpeak[i]){
            statpeak[i] = temp_peak.at(j);
          }
          selected_peaks.push_back(temp_peak.at(j));
          selected_time.push_back(timeg.at(j));
          thismanypointsHigh++;
          totalmany++;
          if(temp_peak[j]<lowerThreshold){
            strikes++;
          }
        }

        if(statpeak[i]>too_big){
          waitingInterval = peakPosition[i]+waiting;
          discard_this[i] = true;
        }
        if(strikes>=maxHits){
          //         cout << "strikes " << strikes << ".... " << eventNow << endl;
          discard_this[i] = true;
        }

        // for get_wave_form
        if(peakPosition.at(i) - mean_before>=0 && peakPosition.at(i)+mean_after<memorydepth*dtime){
          for(Int_t j = peakPosition.at(i)/dtime - mean_before/dtime; j <= peakPosition.at(i)/dtime+mean_after/dtime; j++){
            // temp_waveforms[i].push_back(temp_peak.at(j));
            temp_waveforms[i].push_back(ch.wvf[j]);
            //           if(temp_waveforms[i].size()>200/dtime && temp_waveforms[i].size()<400/dtime){
            //             if(temp_peak.at(j)<-20) notAGoodWaveform[i]=true;
            //           }
          }
        }
        else{
          notAGoodWaveform[i] = true;
        }

        trackHigh[i] = thismanypointsHigh;
        trackLow[i] = thismanypointsLow;




        if(snap()){

          cout << "charge = " << statcharge[i]*dtime << " at " << peakPosition.at(i) << " with " << thismanypointsLow << " + " << thismanypointsHigh << " or this " << thismanypointsBase << " discard ? " << discard_this[i] << endl;
        }






      }
      for(Int_t i = 0; i<npeaks; i++){
        if(discard_this[i]==false){

          peakPerEvent++;
          statcharge[i] = dtime*statcharge[i];
          charge = dtime*charge;
          hcharge->Fill(statcharge[i]);
          hnobase->Fill(statcharge[i]);
          hstat->Fill(statcharge[i]);
          treeCharge = statcharge[i];
          treePeak = statpeak[i];
          ptsHigh = trackHigh[i];
          ptsLow = trackLow[i];
          charge_status += " charge = " + to_string(statcharge[i]);
          tout->Fill();
          //             if (charge<-800){
          //               cout << "\n\n------>event = " << eventNow << " " << charge << " " << peakPosition[i] << endl;
          //             }
          //             if(statcharge[i]>10000 && statcharge[i]<11000 && snap()){
          //                 cout << ".............................. charge = " << statcharge[i] << " at " << peakPosition.at(i) << " event = " << eventNow << endl;
          //             }

          if(get_wave_form){
            Double_t newbase = 0;
            if(notAGoodWaveform[i]==false){
              if(statcharge[i]>=delta/deltaminus && statcharge[i]<=delta*deltaplus){
                naverages++;
                valid = true;
                high_cut = mean_before+mean_after;
                for(Int_t j = low_cut/dtime; j<high_cut/dtime; j++){
                  if(temp_waveforms[i][j]>val_cut){
                    naverages--;
                    valid = false;
                    break;
                  }
                }
              }
              else{
                valid = false;
              }
              for(Int_t j = 0; j<(int)waveforms.size(); j++){
                if(j<(int)temp_waveforms[i].size()){

                  // waveforms[j] = temp_waveforms[i][j];
                  waveforms[j] = temp_waveforms[i][j];
                  sample.wvf[j]=waveforms[j];
                  if(valid) mean_waveforms[j]+=waveforms[j];
                }
                else{
                  waveforms[j] = 0;
                  sample.wvf[j]=0;
                  mean_waveforms[j]+=waveforms[j];
                }

              }
              wvfcharge = statcharge[i];
              sample.event = eventNow;
              sample.charge = statcharge[i];
              sample.selection = valid;
              twvf->Fill();

              //gwaveforms[i] = new TGraph(waveforms.size(),&timeg[0],&waveforms[0]);
              //gm->Add(gwaveforms[i],"LP");
            }
          }


        }
        else{
          //             charge_status += " Discarded - > " + to_string(statcharge[i]*dtime) + " ";
          //             charge_status += " strikes - > " + to_string(trackStrikes[i]) + " ";
        }
      }






    }

    // ____________________________________________________________________________________________________ //
    void searchForPeaks(){
      // For each peak found, we i am looking for the maximum above tolerance
      Int_t n = (int)peak_smooth.size();
      Int_t npeaks = 0;

      higherValue = 0;
      tempPosition = 0;
      Bool_t rebase = false;
      Bool_t normalbase = true; // make it true to calculate
      Double_t gapstart = baselineTime;
      Double_t gapend = 6000;
      Double_t actual_finish = finish;
      //     for(Int_t i = gapstart/dtime; i<=gapend/dtime; i++){
      //       if(temp_peak[i]>200){
      //         actual_finish = gapstart;
      //       }
      //     }

      threshold = tolerance*stddev;
      if(fix_threshold) threshold = tolerance;
      if(check_selection && ch.selection!=0) return;

      for(Int_t i = 0; i<n; i++){

        if(i>timeLow/dtime && i<(n-timeHigh/dtime) && (i>=start/dtime && i<=actual_finish/dtime)){ //searching out of the beginning
          // if(i>baselineTime/dtime && temp_peak[i]>300){ // not sure why this was here
          // if(i>baselineTime/dtime){
          //   break;
          // }
          if(peak_smooth[i]>(mean+threshold) && peak_smooth[i-1]<(mean+threshold) && (i<=gapstart || i>=gapend)){
            npeaks++;

            peakMax.push_back(peak_smooth.at(i));
            peakPosition.push_back(timeg.at(i));

            if(snap()){
              cout << "npeaks = " << npeaks << " " << peakPosition.at(npeaks-1) << " " << peakMax.at(npeaks-1) << " " << eventNow << endl;

            }
            //           if(i<=baselineTime/dtime){
            //               rebase = true;
            //           }
            //           if(rebase && normalbase){
            //             rebase =false;
            //             normalbase=false;
            //             Int_t binmax = hbase->GetMaximumBin();
            //             Double_t newB1 = hbase->GetXaxis()->GetBinCenter(binmax);
            //             Double_t newB2 = hbase->GetMean();
            //             Double_t newB = (newB1<newB2)? newB1 : newB2;
            //
            //             for(Int_t k = 0; k<temp_peak.size();k++){
            //               temp_peak[k] = temp_peak[k]-newB;
            // //               mean = newB;
            //               if(eventNow==my_events[nshow-1]){
            // //                 cout << eventNow << " " << k << " " << temp_peak[k] << " " << newB << " " << hbase->GetMean() << endl;
            //               }
            //             }
            //           }

          }



        }
        else if(i>finish/dtime){
          break;
        }
      }

    }

    vector<Double_t> delay_line(vector<Double_t> v, Double_t delay_time){
      if(delay_time==0) return v;
      vector<Double_t> res(v.size());
      for(int i=0; i<(int)v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
      }
      return res;
    }

    // ____________________________________________________________________________________________________ //
    void lookWaveform(){

      peakPerEvent = 0;

      vector<Double_t> ma_to_shift = movingAverage(&temp_peak[0],pre_filter);
      vector<Double_t> shifted=delay_line(ma_to_shift, shifter);//cusp(MIN[i], h);
      smoothWithMovingAvarage(shifted); // first calculate avarage


      searchForPeaks(); //search for peaks with the moving avarage
      if(snap()){ // if the events match, create the graphs

        g_smooth = new TGraph(timeg.size(),&timeg[0],&peak_smooth[0]);
        g_normal = new TGraph(timeg.size(),&timeg[0],&temp_peak[0]);
        //         cout << mean << " " << stddev << endl;

      }
      integrateSignal();

      if(snap()){
        drawMySamples();
        cout << "Event " << eventNow << " total of peaks: " << peakPerEvent << endl;
      } // draw sample graphs
    }

    // ____________________________________________________________________________________________________ //
    void makeHistogram(string filename){


      if(matching == true){
        getMySample();

        if(shift==0){
          shift = (to-from)/2; //make sure to not use "2." we want an integer here
        }
        if(area_off==0){area_off = to - from;}

        cout << area_off << " " << shift << endl;
      }


      if(creation){
        tout->Branch("value",&value,"value/D");
        tout->Branch("desv",&desv,"desv/D");
        tout->Branch("channel",&channel,"channel/D");
        tout->Branch("treeCharge",&treeCharge,"treeCharge/D");
        tout->Branch("treePeak",&treePeak,"treePeak/D");
        tout->Branch("ptsLow",&ptsLow,"ptsLow/D");
        tout->Branch("ptsHigh",&ptsHigh,"ptsHigh/D");
        tout->Branch("strikes",&strikes,"strikes/D");
        //         creation = false;

        twvf->Branch(Form("Ch%i",channel),&sample,sample.tobranch.c_str());
      }

      cout << "reading: " << filename << endl;
      string rootfile = filename + ".root";

      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");

      TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
      bch->SetAddress(&ch);
      Int_t nentries = t1->GetEntries();

      TH1D *hdark = new TH1D("hdark","",1000,-10000,30000);

      //_____________________________ Start creating random data to check _____________________________ //
      TRandom *rmd = new TRandom();

      my_events[0] = 0;
      for(Int_t i = 1; i<nshow; i++){
        //         my_events[i] = static_cast<Int_t>(rmd->Uniform(my_events[i-1]+1,my_events[i-1]+50));
        my_events[i] = i;
      }
      sort(my_events,my_events+nshow);



      eventNow = 0;
      aux_events = 0;
      //_____________________________ End creating random data to check _____________________________ //



      // ____________________ Start resetting globals ____________________ //
      mean = 0;
      stddev = 0;
      midpoint = 0; // midpoint that takes the avarage
      width = 0;
      sum=0;
      charge = 0;
      discard = false;
      eventNow = 0;
      aux_events = 0;
      teste = 0;
      fillg = true;
      hsample->Reset();
      counter = 0;
      // ____________________ Finish resetting globals ____________________ //




      Double_t aux = 0;
      Int_t j = 0;

      temp_peak.clear();
      temp_peak.resize(memorydepth);
      timeg.clear();
      peak_smooth.clear();
      peakPosition.clear();
      peakMax.clear();
      selected_peaks.clear();
      selected_time.clear();


      matched.clear();

      ymean.resize(nsample);
      ynormal.resize(nsample);
      xmean.resize(nsample);
      for(Int_t i = 0; i<(int)ymean.size(); i++){
        ymean.at(i) = 0;
        xmean.at(i) = 0;
      }


      hbase->Reset();
      hstat->Reset();
      hbase_smooth->Reset();
      hcharge->Reset();
      hzero->Reset();
      hnobase->Reset();

      mean_waveforms.clear();
      mean_waveforms.resize(memorydepth,0);
      naverages = 0;

      for(Int_t i = 0; i<memorydepth; i++){
        timeg.push_back(dtime*i);
      }

      if(just_a_test){nentries = just_this;}
      //     aux = 0;

      for(Int_t i = 0; i<nentries; i++){
        bch->GetEvent(i);
        DENOISE dn;
        if(filter>0)dn.TV1D_denoise<Double_t>(&ch.wvf[0],&temp_peak[0],memorydepth,filter);
        else{
          for(Int_t i = 0; i<memorydepth; i++){
            temp_peak[i] = ch.wvf[i];
          }
        }
        // for(Int_t i = 0; i<memorydepth; i++){
        //   if((i<=5000/dtime) && abs(temp_peak[i])<6){
        //     hbase->Fill(temp_peak[i]);
        //   }
        // }


        aux = ch.event;
        if(static_cast<Int_t>(eventNow)%200==0){
          cout << eventNow << "\r" << flush;
        }
        eventNow =  ch.event;
        // here is were all the calculation is done !!!
        lookWaveform();

        if(snap() && aux_events+1<nshow){
          aux_events = aux_events+1;
        }

        hbase->Reset();
        hbase_smooth->Reset();

        charge_status = "";
        temp_peak.clear();
        temp_peak.resize(memorydepth);
        peak_smooth.clear();
        peakPosition.clear();
        peakMax.clear();
        selected_peaks.clear();
        selected_time.clear();
        baseCharge = 0;
        matched.clear();

      }


      cout << "____________________________ " << counter/nsample << " \n \n \n" << endl;

      if(get_wave_form){
        for(Int_t i = 0; i<memorydepth; i++){
          //             cout << i << " " << mean_waveforms[i] << " ";

          mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
          //             cout << mean_waveforms[i] << endl;

        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;

      }

      TCanvas *cwvf = new TCanvas();
      cwvf->cd();
      //     gm->Draw("A");

      TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);

      if(get_wave_form){
        //         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
      }

      fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

      f1->Close();

    }

    // ____________________________________________________________________________________________________ //
    void smoothWithMovingAvarage(vector<Double_t> shifted){

      Int_t n = shifted.size();
      if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
      }
      else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
      }


      for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-midpoint) || interactions == 0){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          peak_smooth.push_back(shifted.at(i));
        }
        else if(i>cut/dtime){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+shifted.at(j);
            //                 cout << sum << endl;
          }
          peak_smooth.push_back(sum/width);

          if(timeg.at(i)<=baselineTime && abs(peak_smooth.at(i))<baseLimit){
            hbase_smooth->Fill(peak_smooth.at(i));
          }

        }
        else{
          peak_smooth.push_back(0);
        }



        sum=0;
      }
      mean = hbase_smooth->GetMean();
      stddev = hbase_smooth->GetStdDev();



    }
    // ____________________________________________________________________________________________________ //
    void drawMySamples(){

      string sampleName = "ev_" + to_string(my_events[aux_events])+"_"+to_string(channel);

      // TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str(),1920,0,700,500);
      // this is not working when saving
      TCanvas *c1 = new TCanvas();

      c1->cd(1);
      g_smooth->SetLineColor(kRed);
      g_smooth->SetLineWidth(3);

      g_normal->SetLineColor(kBlue);
      g_normal->SetTitle(charge_status.c_str());

      g_normal->Draw("AL");
      g_smooth->Draw("L SAME");

      TLine *lmean = new TLine(timeLimit,mean,memorydepth*dtime,mean);
      TLine *ldev = new TLine(timeLimit,mean+threshold,memorydepth*dtime,mean+threshold);

      lmean->SetLineColor(kGreen);
      ldev->SetLineColor(kGreen);

      lmean->SetLineWidth(2);
      ldev->SetLineWidth(2);

      lmean->Draw();
      ldev->Draw();
      Int_t n = selected_peaks.size();
      if(n!=0){
        g_points = new TGraph(n,&selected_time[0],&selected_peaks[0]);
        g_points->SetMarkerColor(kBlack);
        g_points->Draw("P* SAME");
      }




      // Some fancy drawing now;
      if(matchingAreas.size()!=0){
        Int_t nAreas = matchingAreas.size(); //always two points per area of course
        TLine *larea;
        Double_t max = g_normal->GetYaxis()->GetXmax();
        for(Int_t i = 0; i<nAreas; i++){
          larea = new TLine(matchingAreas[i],0,matchingAreas[i],max);
          cout << "\t\t\t" << matchingAreas[i] << " ";
          if((i+1)%2==0)cout << "\n" << endl;
          larea->SetLineColor(kRed);
          larea->SetLineWidth(2);
          larea->Draw();
        }
      }

      fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");


    }



    // ____________________________________________________________________________________________________ //
    Bool_t snap(){
      if(eventNow == my_events[aux_events]){
        return true;
      }
      else{
        return false;
      }
      return false; // just to be sure oO
    }

    void getMySample(){
      TFile *fsample = new TFile("mysample.root","READ");
      if (fsample->IsZombie()) {
        matching = false;
        std::cout << "\n\n\n\n\n Error opening file \n\n\n\n\n" << std::endl;
        return;
      }
      TTree *tsample = (TTree*)fsample->Get("t1");

      Double_t values;
      tsample->SetBranchAddress("values",&values);
      Int_t nsample = tsample->GetEntries();
      for(Int_t i = 0; i < nsample; i++){
        tsample->GetEntry(i);
        mySample.push_back(values);
      }
      fsample->Close();
      return;

    }

    Bool_t checkAreas(Double_t totalmany){ //return true if the points did not touched the area...
      Int_t nAreas = matchingAreas.size(); //always two points per area of course
      Int_t n = selected_time.size();
      for(Int_t i = n; i>n-totalmany; i--){
        for(Int_t j = 0; j<nAreas; j++){
          if(selected_time[i-1]>=matchingAreas[j] && selected_time[i-1]<=matchingAreas[j+1]){
            return true;
          }
          j++;
        }
      }
      return false;
    }

    void makeSimpleHistogram(string filename){


      sphe_charge = sphe_charge_ch0; // wave0
      sphe_charge2 = sphe_charge2_ch0; // wave0
      delta = sphe_charge2 - sphe_charge;
      sphe_std = sphe_std_ch0;

      if(creation){
        twvf->Branch(Form("Ch%i",channel),&sample,sample.tobranch.c_str());
      }
      cout << "reading: " << filename << endl;
      string rootfile = filename + ".root";

      TFile *f1 = new TFile(rootfile.c_str(),"READ");
      TTree *t1 = (TTree*)f1->Get("t1");

      TBranch *bch = t1->GetBranch(Form("Ch%i",channel));
      bch->SetAddress(&ch);
      Int_t nentries = t1->GetEntries();
      Double_t charge = 0;
      Bool_t noise = false;
      Int_t noise_hits = 0;
      Double_t max = -1e12;
      mean_waveforms.clear();
      timeg.clear();
      mean_waveforms.resize(memorydepth,0);
      naverages = 0;


      // ofstream ftmp;
      // ftmp.open("valid_events.log",ios::out);
      if(just_a_test){nentries = just_this;}
      for(Int_t i = 0; i<nentries; i++){
        bch->GetEvent(i);
        DENOISE dn;
        dn.TV1D_denoise<Double_t>(&ch.wvf[0],&temp_peak[0],memorydepth,filter);

        noise = false;
        max = -1e12;
        for(Int_t j = start/dtime; j<finish/dtime; j++){
          if(withfilter) charge += temp_peak.at(j);
          else charge += ch.wvf[j];
          if(temp_peak[j]>=max){
            max = temp_peak[j];
          }
          if(temp_peak[j]>too_big) noise=true;
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }
        for(Int_t j = 0; j<start/dtime-800/4; j++){
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }
        for(Int_t j = finish/dtime+800/4; j<memorydepth; j++){
          if(temp_peak[j]<lowerThreshold){
            noise_hits++;
            if(noise_hits>=maxHits){
              noise = true;
              break;
            }
          }
        }



        if(noise==false){
          if(check_selection && ch.selection!=0){
          }
          else{
            valid = false;
            for(Int_t j = 0; j<memorydepth; j++){
              sample.wvf[j] = ch.wvf[j];
            }
            if(charge*dtime>=delta/deltaminus  && charge*dtime<=delta*deltaplus){
              valid = true;
              for(Int_t j = 0; j<memorydepth; j++){
                mean_waveforms[j]+=ch.wvf[j];
              }
              naverages++;
            }

            wvfcharge = charge*dtime;
            sample.charge = wvfcharge;
            sample.selection = valid;
            twvf->Fill();
            hcharge->Fill(charge*dtime);
            // cout << charge*dtime << endl;
            // ftmp << i << "\n";
          }
        }
        charge=0;
      }
      // ftmp.close();

      if(get_wave_form){
        for(Int_t i = 0; i<memorydepth; i++){
          timeg.push_back(dtime*i);
          mean_waveforms[i] = mean_waveforms[i]/(naverages == 0 ? 1 : naverages);
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;

      }


      TGraph *gmean = new TGraph(timeg.size(),&timeg[0],&mean_waveforms[0]);
      if(get_wave_form){
        //         fwvf->WriteObject(cwvf,Form("all valid waveforms of ch%i",channel),"TObject::kOverwrite");
        fwvf->WriteObject(gmean,Form("mean_ch%i",channel),"TObject::kOverwrite");
      }
      fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),channel),"TObject::kOverwrite");

      f1->Close();

    }

    template<class T>
    vector<Double_t> movingAverage(T* v, Int_t myinte){

      Int_t n = memorydepth;
      vector<Double_t> res(n,0);
      if(myinte==0) {
        for (Int_t i = 0; i < n; i++) {
          res[i] = v[i];
        }
        return res;
      }
      if(myinte%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = myinte/2+1;    //midpoint will be 5 here
        width = myinte+1;
      }
      else{
        midpoint = (myinte-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = myinte;
      }


      for(Int_t i = 0; i < n; i++){

        if(i<midpoint || i>(n-midpoint)){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          res[i] = v[i];
        }
        else if(i>cut/dtime){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];
            //                 cout << sum << endl;
          }
          res[i] = (sum/width);


        }
        else{
          v[i] = (0);
        }

        sum=0;
      }
      return res;

    }



};
