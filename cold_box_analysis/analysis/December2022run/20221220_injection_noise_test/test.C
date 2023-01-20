TTree *t1 = new TTree("t1","t1");
TEventList *lev;
void test(){

  Double_t x = 0;
  Double_t y = 0;
  t1->Branch("x",&x);
  t1->Branch("y",&y);
  for(Int_t i = 0; i < 500; i++){
    x = i/10.;
    y = x*x;
    t1->Fill();
  }

  t1->Draw(">>lev","Entry$>100");
  lev = (TEventList*)gDirectory->Get("lev");

  TEventList *ltemp = new TEventList("ltemp","ltemp");
  ltemp->Enter(105);
  lev->Subtract(ltemp);
  // lev->Reset();
  t1->Draw(">>lev","");

}
