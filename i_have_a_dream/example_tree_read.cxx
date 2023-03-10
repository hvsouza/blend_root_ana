//
// Usage: root -b -q example_tree_read.cxx
//        root -b -q -e '.L MY_DATA.hxx' example_tree.root -e 't1->Print();'
//

#include "../cold_box_analysis/class/MYCODES.h"

#include "TFile.h"
#include "TTree.h"

#include <iostream>

void example_tree_read() {
  MY_DATA *ch0 = 0, *ch1 = 0; // note: they must be initialized
  TFile *f1 = TFile::Open("example_tree.root", "READ");
  if ((!f1) || f1->IsZombie()) { delete f1; return; } // just a precaution
  TTree *t1; f1->GetObject("t1", t1);
  if (!t1) { delete f1; return; } // just a precaution
  // std::cout << "\n"; t1->Print();
  // t1->SetBranchAddress("ch0.", &ch0);  // note: the name's last character is a "."
  // t1->SetBranchAddress("ch1.", &ch1);  // note: the name's last character is a "."
  TBranch *b1 = t1->GetBranch("ch0.");
  TBranch *b2 = t1->GetBranch("ch1.");
  b1->SetAddress(&ch0);
  b2->SetAddress(&ch1);
  for(Long64_t i = 0; i < t1->GetEntries(); i++) {
    b1->GetEntry(i);
    b2->GetEntry(i);
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
  delete f1; // automatically deletes the "t1", too
  delete ch0;
  delete ch1;
}
