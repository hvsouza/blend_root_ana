//
// Usage: root -b -q example_tree_read.cxx
//        root -b -q -e '.L MY_DATA.hxx' example_tree.root -e 't1->Print();'
//

#include "MY_DATA.hxx"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

void example_tree_read() {
  MY_DATA *ch0 = 0, *ch1 = 0; // note: they must be initialized
  TFile *fread = TFile::Open("example_tree.root", "READ");
  if ((!fread) || fread->IsZombie()) { delete fread; return; } // just a precaution
  TTree *t1; fread->GetObject("t1", t1);
  if (!t1) { delete fread; return; } // just a precaution
  // std::cout << "\n"; t1->Print();
  t1->SetBranchAddress("ch0.", &ch0);  // note: the name's last character is a "."
  t1->SetBranchAddress("ch1.", &ch1);  // note: the name's last character is a "."
  for(Long64_t i = 0; i < t1->GetEntries(); i++) {
    t1->GetEntry(i);
    if (!(ch0 && ch1)) continue; // just a precaution
    std::cout << "\nEntry : " << i << "\n";
    // ch0
    std::cout << "ch0 : " << ch0->peak << " " << ch0->selection << " " << ch0->npts << " :";
    for(Int_t j = 0; j < ch0->npts; j++) std::cout << " " << ch0->wvf[j];
    std::cout << "\n";
    // ch1
    std::cout << "ch1 : " << ch1->peak << " " << ch1->selection << " " << ch1->npts << " :";
    for(Int_t j = 0; j < ch1->npts; j++) std::cout << " " << ch1->wvf[j];
    std::cout << "\n";
  }
  std::cout << "\n";
  // cleanup
  delete fread; // automatically deletes the "t1", too
  delete ch0;
  delete ch1;
}
