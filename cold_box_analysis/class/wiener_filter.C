#include "MYCODES.h"


class WIENER{
public:
  // based https://dspcookbook.github.io/optimal-filtering/wiener-filter-2/#3-solution
  // creates filter from source s and reference d (such as noise)
  vector<Double_t> create_filter(vector<Double_t> s, vector<Double_t> d,Int_t filter_size){
    Int_t n = s.size();
    vector<Double_t> w(filter_size);
    vector<Double_t> rxx(filter_size);
    vector<Double_t> rxd(filter_size);
    TH1D *hnoise = new TH1D("hnoise","hnoise",500,0,0);
    for(Int_t i = 0; i<n; i++){
      hnoise->Fill(d[i]);
    }
    Double_t rvv = hnoise->GetStdDev();
    rvv = rvv*rvv;
    // evaluate rxx (and rxd ?)
    for(Int_t i = 0; i<filter_size; i++){

      for(Int_t j = 0; j<n-i; j++){
        rxx[i] += s[j]*s[j+i];
        rxd[i] += s[j]*d[j+i];
      }
      rxx[i]=rxx[i]/(n-i);
      // rxd[i]=rxd[i]/(n-i);
      if(i==0)rxd[i]=rxx[i]-rvv;
      else rxd[i] = rxx[i];
    }

    TMatrixD Rxx(filter_size,filter_size);
    Int_t aux=0;
    Int_t iteractor;
    for(Int_t i=0; i<filter_size; i++){
      iteractor = aux;
      for(Int_t j = 0; j<filter_size; j++){
        if(iteractor==filter_size) iteractor = 0;
        Rxx[i][j] = rxx[iteractor];
        // cout << Rxx[i][j] << " ";
        iteractor++;
      }
      // cout << "\n";
      aux--;
      if(aux==-1) aux = filter_size-1;
        
    }
    cout << "Inverting Rxx matrix... " << endl; Rxx.Invert(); 
    vector<Double_t> x(filter_size);
    for(Int_t i = 0; i<filter_size; i++){
      x[i] = i;
      for(Int_t j = 0; j<filter_size; j++){
        w[i] += Rxx[i][j]*rxd[j];
      }
      // cout << w[i] << endl;
    }
    
    TCanvas *c2 = new TCanvas();
    TGraph *gtest = new TGraph(filter_size,&x[0],&rxd[0]);
    gtest->Draw("ALP");

    // hnoise->Draw();
    
    return w;
  }

  vector<Double_t> filter(vector<Double_t> s, vector<Double_t> w){
    Int_t n = s.size();
    Int_t m = w.size();
    vector<Double_t> res(n);
    for(Int_t i = 0; i<n; i++){
      // for(Int_t j = i-m/2; j<i+m/2; j++){
      //   if(j<0) continue;
      //   if(j>n) break;
        
      //   res[i]+=s[j]*w[m/2-(j-i)];
      //   // cout << m/2-(j-i) << endl;
      // }

      for(Int_t j = 0; j<m; j++){
        res[i] += s[i]*w[j];
      }
      // cout << res[i] << endl; 
    }
    return res;
  }
  
};
