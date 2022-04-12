// This is an example base on tree2.C tutorial code

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TMath.h"

// this needs to be constant (for writting or reading)
// the value is already set to the maximum of this tree_example.root file npoints variable
// my problem is: suppose I have a second root file in which npoints goes up to 35, how can I change this easily? 
const Int_t MAXMEC = 25;

class DATA{
public:
  Int_t npoints;
  Float_t val[MAXMEC];
  Float_t x[MAXMEC];
};


void tree2w()
{
   //create a Tree file tree_example.root
   TFile f("tree_example.root","recreate");
   TTree t2("t2","t2");
   DATA data;
   t2.Branch("npoints",&data.npoints);
   t2.Branch("myval",data.val,"myval[npoints]/F");
   t2.Branch("x",data.x,"x[npoints]/F");

   // creating a first event
   data.npoints = MAXMEC; // this works if changed to 20;
   for (Int_t i=0;i<data.npoints;i++) {
     data.val[i] = i*i;
     data.x[i] = i/10.;
   }
   t2.Fill();
   // creating a second event   
   data.npoints = 8;
   for (Int_t i=0;i<data.npoints;i++) {
     data.val[i] = i*i*i;
     data.x[i] = i/10.;
   }
   t2.Fill();

   
   t2.Write();
}

void tree_create() {
   tree2w();

}
