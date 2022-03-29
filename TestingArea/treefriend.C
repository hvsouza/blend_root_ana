/// \file
/// \ingroup tutorial_tree
/// \notebook
/// Illustrates how to use Tree friends:
///   - create a simple TTree
///   - Copy a subset of this TTree to a new TTree
///   - Create a Tree Index
///   - Make a friend TTree
///   - compare two TTrees
///   - Draw a variable from the first tree versus a variable
///     in the friend Tree
///
/// You can run this tutorial with:
/// ~~~
///  root > .x treefriend.C  (interpreted via Cling)
///  root > .x treefriend.C+ (executed via ACLIC & the native compiler)
/// ~~~
/// or, variants like:
/// ~~~
///  root > .L treefriend.C+
///  root > CreateParentTree();
///  root > CreateFriendTree();
///  root > CompareTrees();
///  root > DrawFriend();
/// ~~~
///
/// \macro_output
/// \macro_image
/// \macro_code
///
/// \author Rene Brun


#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TTree.h"

class TEST{
public:
Int_t Run, Event;
Float_t x,y,z;
};

void CreateParentTree() {
   // create a simple TTree with 5 branches
   // Two branches ("Run" and "Event") will be used to index the Tree
   TFile *f = new TFile("treeparent.root","recreate");
   TTree *T = new TTree("T","test friend trees");
   TEST test;
   T->Branch("test",&test,"Run/I:Event/I:x/F:y/F:z/F");

   TRandom *r = new TRandom(2);
   for (Int_t i=0;i<10000;i++) {
      if (i < 5000) test.Run = 1;
      else          test.Run = 2;
      test.Event = i;
      test.x = r->Gaus(10,1);
      test.y = r->Gaus(20,2);
      test.z = r->Landau(2,1);
      T->Fill();
   }
   T->Print();
   T->Write();
   delete f;
}
void CreateFriendTree() {
   // Open the file created by CreateParentTree
   // Copy a subset of the TTree into a new TTree
   //   (see also tutorials copytree.C, copytree2.C and copytree3.C)
   // Create an index on the new TTree ("Run","Event")
   // Write the new TTree (including its index)

  TFile *ff = new TFile("treefriend.root","recreate");
  TTree *TF = new TTree("T","test friend trees");
  TEST test;
  TF->Branch("test",&test,"Run/I:Event/I:x/F:y/F:z/F");

   TRandom *r = new TRandom(1);
   for (Int_t i=0;i<10000;i++) {
      if (i < 5000) test.Run = 1;
      else          test.Run = 2;
      test.Event = i;
      test.x = r->Gaus(20,3);
      test.y = r->Gaus(20,2);
      test.z = r->Landau(2,1);
      TF->Fill();
   }
   TF->Print();
   TF->Write();
   delete ff;
}

void CompareTrees() {
   // The two TTrees created above are compared.
   // The subset of entries in the small TTree must be identical
   // to the entries in the original TTree.

   TFile *f = new TFile("treeparent.root");
   TTree *T  = (TTree*)f->Get("T");
   TFile *ff = new TFile("treefriend.root");
   TTree *TF = (TTree*)ff->Get("T");
   TEST test;
   T->SetBranchAddress("test",&test);
   // T->SetBranchAddress("Event",&Event);

   T->AddFriend(TF);
   Int_t fRun,fEvent;
   Float_t fx,fy,fz;

   Long64_t nentries = T->GetEntries();
   Int_t nok = 0;
   for (Long64_t i=0;i<nentries;i++) {
      T->GetEntry(i);
      if (fRun == test.Run && fEvent==test.Event && test.x==fx && test.y==fy &&test.z==fz) {
         nok++;
      } else {
         if (TF->GetEntryWithIndex(test.Run,test.Event) > 0) {
            if (i <100) printf("i=%lld, Run=%d, Event=%d, x=%g, y=%g, z=%g,  : fRun=%d, fEvent=%d, fx=%g, fy=%g, fz=%g\n",i,test.Run,test.Event,test.x,test.y,test.z,fRun,fEvent,fx,fy,fz);
         }
      }
   }
   printf("nok = %d, fentries=%lld\n",nok,TF->GetEntries());

   delete f;
   delete ff;
}

void DrawFriend() {
  // Draw a scatter plot of variable x in the parent TTree versus
  // the same variable in the subtree.
  // This should produce points along a straight line.

   TFile *f  = TFile::Open("treeparent.root");
   TTree *T  = (TTree*)f->Get("T");
   TFile *f2 = TFile::Open("treefriend.root");
   TTree *t2 = (TTree*)f2->Get("T");
   T->AddFriend(t2,"TF",f2);
   T->Draw("test.x:TF.test.x");
}

void treefriend() {
   CreateParentTree();
   CreateFriendTree();
   // CompareTrees();
   DrawFriend();
}
