//
// Usage: root -b -q example_tree_write.cxx
//        root -b -q -e '.L MY_DATA.hxx' example_tree.root -e 't1->Print();'
//

#include "MY_DATA.hxx"

#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

void example_tree_write() {
  Int_t npts, max_npts = 10;
  Int_t nchannels = 2;
  vector<MY_DATA*> ch(nchannels);
  TFile *fnew = TFile::Open("example_tree.root", "RECREATE");
  TTree *t1 = new TTree("t1", "t1");
  for(Int_t i = 0; i < nchannels; i++) t1->Branch(Form("ch%d.",i), &ch[i]); // note: the name's last character is a "."
  for(Int_t i = 0; i < 5; i++) { // 5 events only
    // ch[0]
    npts = 10;//Int_t((max_npts) * gRandom->Rndm() + 0.5); // variable size
    ch[0]->Set_npts(npts);
    ch[0]->peak = i;
    ch[0]->selection = (max_npts) * npts + i;
    for(Int_t j = 0; j < npts; j++) ch[0]->wvf[j] = (max_npts) * npts + j;
    // ch[1]
    npts = Int_t((max_npts) * gRandom->Rndm() + 0.5); // variable size
    ch[1]->Set_npts(npts);
    ch[1]->peak = i;
    ch[1]->selection = (max_npts) * npts + i;
    for(Int_t j = 0; j < npts; j++) ch[1]->wvf[j] = (max_npts) * npts + j;
    //
    t1->Fill();
  }
  t1->Write();
  delete fnew; // automatically deletes the "t1", too
}
