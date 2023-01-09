// ________________________________________ //
// Author: Henrique Souza
// Filename: readingCodes.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"


class DENOISE{
  
public:
    template<class T>
    void TV1D_denoise(T* input, T* output, unsigned int width, const Double_t lambda) {

        unsigned int* indstart_low = (unsigned int*) malloc(sizeof *indstart_low * width);
        unsigned int* indstart_up = (unsigned int*) malloc(sizeof *indstart_up * width);
        unsigned int j_low = 0, j_up = 0, jseg = 0, indjseg = 0, i=1, indjseg2, ind;
        double output_low_first = input[0]-lambda;
        double output_low_curr = output_low_first;
        double output_up_first = input[0]+lambda;
        double output_up_curr = output_up_first;
        const double twolambda=2.0*lambda;
        if (width==1) {output[0]=input[0]; return;}
        indstart_low[0] = 0;
        indstart_up[0] = 0;
        width--;
        for (; i<width; i++) {
            if (input[i]>=output_low_curr) {
                if (input[i]<=output_up_curr) {
                    output_up_curr+=(input[i]-output_up_curr)/(i-indstart_up[j_up]+1);
                    output[indjseg]=output_up_first;
                    while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
                        output_up_curr+=(output[ind]-output_up_curr)*
                            ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
                    if (j_up==jseg) {
                        while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
                            indjseg2=indstart_low[++jseg];
                            output_up_curr+=(output_up_curr-output_low_first)*
                                ((double)(indjseg2-indjseg)/(i-indjseg2+1));
                            while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                            output_low_first=output[indjseg];
                        }
                        output_up_first=output_up_curr;
                        indstart_up[j_up=jseg]=indjseg;
                    } else output[indstart_up[j_up]]=output_up_curr;
                } else
                    output_up_curr=output[i]=input[indstart_up[++j_up]=i];
                output_low_curr+=(input[i]-output_low_curr)/(i-indstart_low[j_low]+1);
                output[indjseg]=output_low_first;
                while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
                    output_low_curr+=(output[ind]-output_low_curr)*
                            ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
                if (j_low==jseg) {
                    while ((output_low_curr>=output_up_first)&&(jseg<j_up)) {
                        indjseg2=indstart_up[++jseg];
                        output_low_curr+=(output_low_curr-output_up_first)*
                            ((double)(indjseg2-indjseg)/(i-indjseg2+1));
                        while (indjseg<indjseg2) output[indjseg++]=output_up_first;
                        output_up_first=output[indjseg];
                    }
                    if ((indstart_low[j_low=jseg]=indjseg)==i) output_low_first=output_up_first-twolambda;
                    else output_low_first=output_low_curr;
                } else output[indstart_low[j_low]]=output_low_curr;
            } else {
                output_up_curr+=((output_low_curr=output[i]=input[indstart_low[++j_low] = i])-
                    output_up_curr)/(i-indstart_up[j_up]+1);
                output[indjseg]=output_up_first;
                while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
                    output_up_curr+=(output[ind]-output_up_curr)*
                            ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
                if (j_up==jseg) {
                    while ((output_up_curr<=output_low_first)&&(jseg<j_low)) {
                        indjseg2=indstart_low[++jseg];
                        output_up_curr+=(output_up_curr-output_low_first)*
                            ((double)(indjseg2-indjseg)/(i-indjseg2+1));
                        while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                        output_low_first=output[indjseg];
                    }
                    if ((indstart_up[j_up=jseg]=indjseg)==i) output_up_first=output_low_first+twolambda;
                    else output_up_first=output_up_curr;
                } else output[indstart_up[j_up]]=output_up_curr;
            }
        }
        /* here i==width (with value the actual width minus one) */
        if (input[i]+lambda<=output_low_curr) {
            while (jseg<j_low) {
                indjseg2=indstart_low[++jseg];
                while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                output_low_first=output[indjseg];
            }
            while (indjseg<i) output[indjseg++]=output_low_first;
            output[indjseg]=input[i]+lambda;
        } else if (input[i]-lambda>=output_up_curr) {
            while (jseg<j_up) {
                indjseg2=indstart_up[++jseg];
                while (indjseg<indjseg2) output[indjseg++]=output_up_first;
                output_up_first=output[indjseg];
            }
            while (indjseg<i) output[indjseg++]=output_up_first;
            output[indjseg]=input[i]-lambda;
        } else {
            output_low_curr+=(input[i]+lambda-output_low_curr)/(i-indstart_low[j_low]+1);
            output[indjseg]=output_low_first;
            while ((j_low>jseg)&&(output_low_curr>=output[ind=indstart_low[j_low-1]]))
                output_low_curr+=(output[ind]-output_low_curr)*
                            ((double)(indstart_low[j_low--]-ind)/(i-ind+1));
            if (j_low==jseg) {
                if (output_up_first>=output_low_curr)
                    while (indjseg<=i) output[indjseg++]=output_low_curr;
                else {
                    output_up_curr+=(input[i]-lambda-output_up_curr)/(i-indstart_up[j_up]+1);
                    output[indjseg]=output_up_first;
                    while ((j_up>jseg)&&(output_up_curr<=output[ind=indstart_up[j_up-1]]))
                        output_up_curr+=(output[ind]-output_up_curr)*
                            ((double)(indstart_up[j_up--]-ind)/(i-ind+1));
                    while (jseg<j_up) {
                        indjseg2=indstart_up[++jseg];
                        while (indjseg<indjseg2) output[indjseg++]=output_up_first;
                        output_up_first=output[indjseg];
                    }
                    indjseg=indstart_up[j_up];
                    while (indjseg<=i) output[indjseg++]=output_up_curr;
                }
            } else {
                while (jseg<j_low) {
                    indjseg2=indstart_low[++jseg];
                    while (indjseg<indjseg2) output[indjseg++]=output_low_first;
                    output_low_first=output[indjseg];
                }
                indjseg=indstart_low[j_low];
                while (indjseg<=i) output[indjseg++]=output_low_curr;
            }
        }
        free(indstart_low);
        free(indstart_up);
    }




    template<class T>
    void movingAverage(T* v, T* res, Int_t myinte = 10){

      Int_t n = memorydepth;
      Int_t midpoint = 0;
      Int_t width = 0;
      Double_t sum = 0;

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
        else{
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];
            //                 cout << sum << endl;
          }
          res[i] = (sum/width);
        }


        // }

        sum=0;

      }

    }
};










class Headers{

public:

  int EventSize;
  int Board_id;
  int Pattern;
  int Channel;
  int EventCounter;
  int TriggerTimeTag;
};

class Read{
  
public:
  
  

  // Created by Henrique Souza, June 2020
  

  
  Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
  Int_t nbits = 14;
  Int_t basebits = nbits;
  Bool_t isBinary = false;
  Bool_t saveFilter = false;
  
  Double_t startCharge = 3300;
  Double_t maxRange = 5000;
  Double_t fast = 500;
  Double_t slow = 17000;
  Double_t filter = 9;
  Double_t exclusion_baseline = 30;
  Double_t exclusion_window = 500;
  Double_t currentTime = 0;
  Double_t timeCicle = TMath::Power(2,31)-1;
  Double_t timeResolution = 8e-9; // 8 ns for 2 ns step, 16 ns for 4 ns step
  Double_t timestamp = 0;
  Double_t deltastamp = 0;
  Double_t temptime = 0;

  string file_extension = ".txt";
  
  Bool_t OnlyOneEvent = false; // Do you want only one event? Choose it wisely 
  Int_t stopEvent = 2000;
  
  Int_t maxEvents = 100000;
    
  Double_t baselineTime = 10000; // time limit to start looking for baseline
  Double_t chargeTime = 18000; // last time to integrate
  Bool_t noBaseline=false;
  // Bool_t noBaseline=true;
  vector<Int_t> channels = {1,2};
  vector<Double_t> exclusion_baselines = {30};
  Int_t nfiles = 1;

  // THIS did not work properly, sometimes it gets stuck in a bad region. 
  // TH1D *hbase = new TH1D("hbase","finding baseline",TMath::Power(2,nbits),0,0);

  TH1D *htests = new TH1D("htests","htests",1000,0,0);
  TH1D *hbase;
  TF1* fbase = new TF1("fbase","gaus(0)",0,TMath::Power(2,nbits));


  vector<Int_t> navg;
  vector<vector<Double_t>> avg;
  vector<Double_t> time;



  string format_time = "hh:mm:ss";
  string format_date = "dd-mmm-yyyy";

  void get_ch_info(string logfilename){
    ifstream logfile;
    logfile.open(logfilename.c_str(),ios::in);

    if(logfile.good() && logfile.is_open()){ // Ok
      // cout << "Reading file " << files << " ... " << endl;
    }
    else{
      cout << "The file " << logfilename << " did not open!!" << endl;
      return;
    }
    string dataname = "";
    Int_t aux = 0;
    // While for read the entire file
    string wave_ref = "wave";
    int wavenum = 0;
    int temp = -1;
    vector<Int_t> tempch;
    while(!logfile.fail()){
        logfile >> dataname;
        // cout << dataname << endl;
        if(logfile.bad() || logfile.fail()){
          break;
        }
        size_t found  = dataname.find(wave_ref);
        found = found + wave_ref.length();
        wavenum = (int)dataname[found] - '0';
        // cout << wavenum << endl;
        if(wavenum != temp){
          if (aux == 0) temp = wavenum;
          tempch.push_back(wavenum);
        }
        else{
          break;
        }
        aux++;
    }
    channels = tempch;
    Int_t nchannels = channels.size();
    Int_t noriginal = exclusion_baselines.size();
    if (noriginal > nchannels){
      if (noriginal-1 >= channels[nchannels-1]){
        for (Int_t j = 0; j < nchannels; j++) {
          exclusion_baselines[j] =  exclusion_baselines[channels[j]];
        }
        exclusion_baselines.resize(nchannels);
      }
    }
    else if (noriginal < nchannels && noriginal != 1) {
      exclusion_baselines.resize(nchannels);
      for(Int_t j = noriginal; j < exclusion_baselines.size(); j++){
        exclusion_baselines[j] = exclusion_baselines[0];
      }
    }
    // for (Int_t j = 0; j < nchannels; j++) cout << j << " " << channels[j] << " " <<  exclusion_baselines[j] << endl;
    logfile.close();
  }
  
  void adc_read_all_data(Bool_t do_get_ch_info = true){

    if (do_get_ch_info) get_ch_info("files.log");
    readFiles("files.log"); //use it like this
    return;
    
  }
  
  void adc_read(string filename){
    
    readFile(filename);
  }
  
  void readFile(string file){
    if(channels.size()>1){
      channels = {1};
    }
    vector<string> file_ch(channels.size());
    string rootfile;
    file_ch[0] = file;

    rootfile = file + ".root";
    TFile *f1;
    TTree *t1;
    Double_t tEvent = 0;
    vector<ADC_DATA> ch(channels.size());
    vector<TBranch*> bch(channels.size());
    avg.resize(channels.size());
    navg.resize(channels.size(),0);
    for(Int_t i = 0; i<channels.size(); i++){
      avg[i].resize(memorydepth,0);
    }
    
    
    f1 = new TFile(rootfile.c_str(),"RECREATE");
    t1 = new TTree("t1","ADC processed waveform");
    for(Int_t i = 0; i<channels.size(); i++){
      bch[i] = t1->Branch(Form("Ch%i",channels[i]),&ch[i],ch[0].tobranch.c_str());
    }
    
    
    f1->Write();
    
    f1->Close();
  
      // _______________________________________________________________________________________________________
      
    readData(file_ch,rootfile,tEvent);
    
    time.resize(memorydepth);
    for(Int_t i = 0; i<memorydepth; i++){
      time[i] = i*dtime;
    }
    for(Int_t i = 0; i<channels.size(); i++){
      for(Int_t j = 0; j<memorydepth; j++){
        avg[i][j] = avg[i][j]/navg[i];
      }
    }
    
    f1 = new TFile(rootfile.c_str(),"UPDATE");
    vector<TGraph*> gavg(channels.size());
    for(Int_t i = 0; i<channels.size(); i++){
      gavg[i] = new TGraph(memorydepth,&time[0],&avg[i][0]);
      f1->WriteObject(gavg[i],Form("average_ch%i",channels[i]),"TObject::kOverwrite");
    }
//     f1->Close();
    
    
  }
  
  void readFiles(string files){
    
    ifstream filesdata;
    filesdata.open(files.c_str(),ios::in);
    vector<string> file_ch(channels.size());
    string rootfile;
    TFile *f1;
    TTree *t1;
    hbase = new TH1D("hbase","finding baseline",TMath::Power(2,basebits),0,TMath::Power(2,nbits));
    Double_t tEvent = 0;
    vector<ADC_DATA> ch(channels.size());
    vector<TBranch*> bch(channels.size());

    
    Bool_t first_file = true;
    avg.resize(channels.size());
    navg.resize(channels.size(),0);
    for(Int_t i = 0; i<channels.size(); i++){
      avg[i].resize(memorydepth,0);
    }

    if(isBinary) file_extension = ".dat";
    else file_extension = ".txt";
    
    if(filesdata.good() && filesdata.is_open()){ // Ok
      cout << "Reading file " << files << " ... " << endl;
    }
    else{ 
      cout << "The file " << files << " did not open!!" << endl;
      return;
    }
    string temp = "";
    // While for read the entire file
    while(!filesdata.fail()){
      for(Int_t i = 0; i<channels.size(); i++){
        filesdata >> temp;
        file_ch[i] = temp;
      }
      // cout << "C1 = " << file_ch[0] << " C2 = " << file_ch[1] << endl;
      if(filesdata.bad() || filesdata.fail()){
        break;
      }
      
      // ______________________________ Create root files at first files only __________________________________
      if(first_file){
        first_file = false;
        
        rootfile = "analyzed.root";
        
        string erase = "rm " + rootfile;
        system(erase.c_str());
        
        f1 = new TFile(rootfile.c_str(),"RECREATE");
        t1 = new TTree("t1","ADC processed waveform");
        for(Int_t i = 0; i<channels.size(); i++){
          bch[i] = t1->Branch(Form("Ch%i",channels[i]),&ch[i],ch[0].tobranch.c_str());
        }
        
        
        f1->Write();
        
        f1->Close();
      }
      // _______________________________________________________________________________________________________
      
      readData(file_ch,rootfile,tEvent);
      
    }
    time.resize(memorydepth);
    for(Int_t i = 0; i<memorydepth; i++){
      time[i] = i*dtime;
    }
    for(Int_t i = 0; i<channels.size(); i++){
      for(Int_t j = 0; j<memorydepth; j++){
        avg[i][j] = avg[i][j]/navg[i];
      }
    }
    
    f1 = new TFile(rootfile.c_str(),"UPDATE");
    vector<TGraph*> gavg(channels.size());
    for(Int_t i = 0; i<channels.size(); i++){
      gavg[i] = new TGraph(memorydepth,&time[0],&avg[i][0]);
      f1->WriteObject(gavg[i],Form("average_ch%i",channels[i]),"TObject::kOverwrite");
    }
//     f1->Close();

    // htests->Draw();
    
  }
  




  
  // This function will read your data and create a root file with the same name
  void readData(vector<string> file_ch,string rootfile, Double_t &tEvent){
    
    TFile *f1 = new TFile(rootfile.c_str(),"UPDATE");
    TTree *t1 = (TTree*)f1->Get("t1");
    
    
    vector<ADC_DATA> ch(channels.size());
    vector<TBranch*> bch(channels.size());
    vector<string> filename_ch(channels.size());
    TimeREADER myTimer;
    
    DENOISE dn;
    

    if(dtime == 2)  timeResolution = 8e-9; // 8 ns for 2 ns step, 16 ns for 4 ns step
    if(dtime == 4)  timeResolution = 16e-9; // 16 ns for 4 ns step
    
    Int_t numberoflines=0;
    Double_t temp;
    vector<Double_t> raw(memorydepth);
    Double_t filtered[memorydepth];
    
    for(Int_t i = 0; i<channels.size(); i++){
      bch[i] = t1->GetBranch(Form("Ch%i",channels[i]));
      bch[i]->SetAddress(&ch[i]);
      filename_ch[i] = file_ch[i] + file_extension;
    }

    
    vector<ifstream> fin(channels.size());
    for(Int_t i = 0; i<channels.size(); i++){
      if(!isBinary) fin[i].open(filename_ch[i].c_str(),ios::in);
      else          fin[i].open(filename_ch[i].c_str(),ios::in | ios::binary);
      // here i check if the data file is open (or exist), if not, I remove the .root files and close the program. 
      if(fin[i].good() && fin[i].is_open()){ // Ok
        cout << "Reading file " << filename_ch[i] << " ... " << endl;
      }
      else{ 
        cout << "File " << filename_ch[i] << " did not open!!" << endl;
        return;
        
      }
    }
    Bool_t closeMyWhile = false; // keep it false
    // While for read the entire file
    Double_t bl = 0;
    Int_t eventFile = 0;
    string headers;
    Double_t init_time = 0;
    uint32_t valbin = 0;
    Headers headbin;
    int nbytes = 4;
    Int_t headers_npoints = 0;
    Int_t headers_nwvfs = 0;
    Bool_t withTimestamp=true;

    string date, time;
    Double_t stamp;
    vector<Double_t> event_time;
    Double_t starting_time = 0;
    Int_t aux_time = 0;
    event_time.resize(memorydepth);
    for(Int_t i = 0; i<channels.size(); i++){
      if(isBinary==false){
        // LECROYWR104MXi-A/� 49455 Waveform
        // Segments 2000 SegmentSize 2502
        // Segment TrigTime TimeSinceSegment1
        // #1 01-Jan-2002 00:32:41 0                 
        // #2 01-Jan-2002 00:32:41 0.0001    
        // continue;
        getline(fin[i],headers);
        // cout << headers << endl;
        fin[i] >> headers >> headers_nwvfs >> headers >> headers_npoints;
        // cout << headers << endl;
        getline(fin[i],headers); // taking extra \r

        if(headers_npoints!=memorydepth && headers_npoints!=memorydepth-1){
          cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ \n\n correct the memorydepth!!! " << endl;
          cout << "current: " << memorydepth << endl;
          cout << "set to: " << headers_npoints << endl; 
          return;
        }
        //           else{
        //               cout << "@@@@@ \n\n ALLLL GOOOOOD!!!" << endl;
        //               return 1;
        //           }
        //             cout << headers << endl;
        //         }
        
        if(withTimestamp){
          getline(fin[i],headers); // reads the header of time stamp
          for(Int_t ln=0;ln<headers_nwvfs;ln++){
            fin[i] >> headers >> date >> time >> stamp;
            // cout << ln <<  " " << headers_nwvfs << " " <<  headers << " " << date << " " <<  time << endl;
            if(ln==0){
              // cout << date << " " << time << endl;
              starting_time = myTimer.timeRead(date,time,format_date,format_time);
              event_time[ln] = starting_time;
              // printf("Starting time = %1.f\n",starting_time);
            }
            else{
              starting_time+=stamp;
              event_time[ln] = starting_time;
            }
          }
          getline(fin[i],headers); // reads the extra \r
        }
         
        getline(fin[i],headers);
        // cout << headers << endl;
      }
    }

    if(isBinary){
      fin[0].read((char *) &headbin, nbytes*6);
      if((headbin.EventSize-24)/2 != memorydepth){
        cout << "ERROR !!!! INCORRECT LENGTH !!! Reading: " << memorydepth << ", folder: " << (headbin.EventSize-24)/2 << endl;
        return;
      }
    }

    fin[0].clear();
    fin[0].seekg(0);
    while(!fin[0].fail() && closeMyWhile == false){ // We can just look for one file, they shold have the same amount of lines anyway!
      Int_t n_reads = 0;    
    
      for(Int_t i = 0; i<channels.size(); i++){
        if(isBinary==false){ 

          
          for(int j = 0; j < memorydepth; j++)
            {
              if(withTimestamp) fin[i] >> timestamp >> temp;
              else fin[i] >> temp;
              if(fin[i].bad() || fin[i].fail()){
                cout << "going ... ";
                break;
              }
              n_reads++;
              raw[j] = temp;
              ch[i].wvf[j] = temp;
              filtered[j] = temp;
            }
          //           cout << "............................ \n";
          //           getline(fin[i],headers);
        }
        else{
          // for(Int_t ln=0;ln<6;ln++){ // 4 bytes (32 bits) for each head (no text) 
            fin[i].read((char *) &headbin, nbytes*6);
            // printf("%d\n",headbin);
          // }
          timestamp = headbin.TriggerTimeTag;

          //           printf("%.0f\n",timestamp);
          for(int j = 0; j < memorydepth; j++)
            {
              fin[i].read((char *) &valbin, 2);
              //             if(j==0) printf("%d\n -- \n",valbin);

              if(fin[i].bad() || fin[i].fail()){
                break;
              }
              // cout << valbin << endl;
              raw[j] = valbin;
              ch[i].wvf[j] = valbin;
              filtered[j] = valbin;
            }
        }

        if((fin[i].bad() || fin[i].fail()) && n_reads<memorydepth-5){
          // cout << "problems ??" << endl;
          break; // giving a 5 points relaxiation 
        }
        if(i==0){
          if(timestamp<0){
            timestamp = timeCicle+timestamp;
          }
          if(timestamp<temptime){
            deltastamp = timestamp+timeCicle-temptime;
          }
          else{
            deltastamp = timestamp - temptime;
          }

          temptime = timestamp;

          if(init_time!=0 && eventFile<maxEvents){
            currentTime = currentTime+deltastamp*timeResolution;
          }
          else{
            init_time = 1;
          }
          ch[i].time = currentTime;
        }

        // if(!isBinary) ch[i].time = event_time[aux_time];
        // printf("time of event = %11f\n",event_time[aux_time]);
        aux_time++;
        if(filter>0) dn.TV1D_denoise<Double_t>(&ch[i].wvf[0],&filtered[0],memorydepth,filter);
        // if(filter>0) dn.TV1D_denoise<Double_t>(&raw[0],&ch[i].wvf[0],memorydepth,filter);
        if(saveFilter==true){
          for(Int_t l = 0; l<memorydepth; l++){
            ch[i].wvf[l] = filtered[l];
          }
        }

        if(exclusion_baselines.size() != 1){
          exclusion_baseline = exclusion_baselines[i];
        }
        bl = baseline(filtered,ch[i].selection,i,tEvent);
        // bl = baseline(ch[i].wvf,ch[i].selection,i,tEvent);
        // if(bl==-9999) cout << i << " " << tEvent << endl;
        getvalues(i,ch[i],filtered,bl);
        ch[i].event = tEvent;
        
        numberoflines++;
        
        
        
      }
      if((fin[0].bad() || fin[0].fail()) && n_reads<memorydepth-5){
        // cout << "problems ??" << endl;
        break; // giving a 5 points relaxiation 
      }
      // if(fin[0].eof()){
        // numberoflines--;
        // break;
      // }
      if(eventFile<maxEvents){
        t1->Fill();
        tEvent+=1;
      }

      if(OnlyOneEvent == true && eventFile==stopEvent-1){
        closeMyWhile=true;
      }
      eventFile++;
      
    }
    cout << tEvent << " events in " << numberoflines << " lines " << endl;
    
    f1->WriteObject(t1,"t1","TObject::kOverwrite");
    
    f1->Close();
  
  }
  
  

  void getvalues(Int_t &nch,ADC_DATA &ch,Double_t filtered[],Double_t bl){
    
    ch.peak =0;
    Double_t fastcomp = 0;
    Double_t slowcomp = 0;
    ch.charge=0;
    navg[nch]++;
    for(Int_t i = 0; i<memorydepth; i++){
      ch.wvf[i] = ch.wvf[i]-bl;
      filtered[i] = filtered[i]-bl;
      avg[nch][i]+=ch.wvf[i];
//       cout << i << " " << ch.wvf[i] << endl;
      if(i>=startCharge/dtime && i<chargeTime/dtime){
        ch.charge+=filtered[i]*dtime;
        if(i <= maxRange/dtime){
          if(ch.peak==0){ ch.peak = filtered[i]; ch.peakpos = i*dtime;}
          else if(ch.peak<filtered[i]){ch.peak=filtered[i];ch.peakpos = i*dtime;}
        }
        if(i<(startCharge+fast)/dtime){
          fastcomp+=ch.wvf[i];
        }
        if(i<(startCharge+slow)/dtime){
          slowcomp+=ch.wvf[i];
        }

      }
    }
//     cout << fastcomp << " " << slowcomp << endl;
    ch.fprompt = fastcomp/slowcomp;
  }
  
  
  Double_t baseline(Double_t v[],Int_t &selection, Int_t idx, Int_t mevent){
    if(noBaseline) return 0;
    Double_t result = 0;
    hbase->Reset();
    for(Int_t i=0; i<baselineTime/dtime; i++) hbase->Fill(v[i]);
    Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());
    Double_t hmean = hbase->GetMean();
    Double_t hstd = hbase->GetStdDev();
    
    // fbase->SetParameters(1000,res0,hstd);
    // hbase->Print();
    
    // // to debug ! Note, fitting the baseline makes its about 10 times slower..

    Bool_t changed_mean = false;
    // if(hstd/hmean>0.1) hstd = 0.05*hmean; // in the case the standard deviation is bigger then X%, we correct it
    // if(res0>hmean+hstd || res0<hmean-hstd){ // res0 is way too out of the mean, so:

    // TFitResultPtr r = hbase->Fit("fbase","WQ0");
    //   Int_t fitStatus = r;
    //   // if(fitStatus==0){
    //     hmean = fbase->GetParameter(1);
    //     hstd = fbase->GetParameter(2);
    //     res0 = hmean;
    //     changed_mean = true;
    //   // }
    // }
    
    // if(idx == 0 && mevent == 3294){
    //   cout << hmean << " " << hstd << " " << res0 << endl;
    // }
    
    Double_t bins=0;
    for(Int_t i=0; i<baselineTime/dtime;){
      if(v[i] > res0 + exclusion_baseline || v[i]<res0 - exclusion_baseline) {
        i+=exclusion_window/dtime;
      }
      else{
        result += v[i];
        bins+=1;
        i++;
      }
    }
    if(bins>0)result/=bins;
    if(bins > (baselineTime/dtime)/3.){
      selection = 0;
      // // You can use this to debug. If selection == 0, there should not have events here.
      // // this means that res0 should be pretty much the average for a good baseline. 
      // if(res0>hmean+hstd || res0<hmean-hstd){
      //   cout << "PROBLEEEMMM " << endl;
      //   cout << result << endl;
      //     cout << res0 << endl;
      //     cout << hmean << endl;
      //     cout << hstd << endl;
      //     // cout << fitStatus << endl;
      //     result = -9999;
      // }


      return result;
    }
    else{
      // cout << result << " " << res0 << endl;
      // htests->Fill(hstd/hmean);
      // if(hstd/hmean<0.4) cout << idx << " " << mevent << endl;
      if(changed_mean==false)selection = 1;
      else selection = 2;
      return res0;
    }
  }
};






















































