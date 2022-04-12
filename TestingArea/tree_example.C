// This is an example base on tree2.C tutorial code


// this needs to be constant (for writting or reading)
// the value is already set to the maximum of this tree_example.root file npoints variable
// my problem is: suppose I have a second root file in which npoints goes up to 35, how can I change this easily? 

#include <iostream>

using namespace std;
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"

class DATA{
public:
  Int_t npoints = 0;
  Float_t *val = nullptr;
  Float_t *x = nullptr;
  DATA(const Int_t nmax){
    val = new Float_t[nmax];
    x = new Float_t[nmax];
  }
};


void tree2r()
{
  TFile *f = new TFile("tree_example.root","READ");
   TTree *t2 = (TTree*)f->Get("t2");

   DATA *data = new DATA(30);

   t2->SetBranchAddress("npoints",&data->npoints);
   t2->SetBranchAddress("myval",data->val);
   t2->SetBranchAddress("x",data->x);

   
   Long64_t nentries = t2->GetEntries();
   vector<TGraph*> g(nentries);
   for (Long64_t i=0;i<nentries;i++) {
     t2->GetEntry(i);

     g[i] = new TGraph(data->npoints,data->x,data->val);  
   }
   
   TCanvas *c1 = new TCanvas("c1","c1",600,800);
   c1->cd(1);
   g[0]->Draw("ALP");
   g[1]->Draw("SAME LP");
}

void tree_example() {
   tree2r();
}
