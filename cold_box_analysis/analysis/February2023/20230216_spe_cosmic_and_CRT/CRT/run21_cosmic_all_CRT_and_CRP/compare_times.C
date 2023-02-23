#define memorydepth 2500
#include "/home/henrique/Dropbox/APC_Paris/Root/cold_box_analysis/class/MYCODES.h"

void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.15,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

class DATA {
  public:
    vector<Double_t> events;
    vector<Double_t> time;
    Double_t ev;
    Double_t ti;
    Int_t n = 0;
    // DATA();
    // virtual ~DATA();

    void fill(TH2D *h){
      for(auto mt: time){
        h->Fill(mt,0.5);
      }
    }

};

void read_data(DATA &EV, const char* file){
  ifstream f;
  f.open(file,ios::in);

  string header;
  getline(f,header);
  char comma{};
  while(f >> EV.ev >> comma >> EV.ti){
    EV.events.push_back(EV.ev);
    EV.time.push_back(EV.ti);
    EV.n++;
  }
}
void compare_times(){

  DATA PDS;
  DATA CRT;
  read_data(PDS,"PDS_TriggerTimes.txt");
  read_data(CRT,"CRT_TriggerTimes.txt");

  vector<TH2D *> h(3);
  for(Int_t i = 0; i < 3; i++){
    h[i] = new TH2D(Form("h%d",i),"",132*10,0,132,1,0,1);
  }
  PDS.fill(h[0]);
  CRT.fill(h[1]);
  TCanvas *c = new TCanvas("c", "c",1920,0,1920,1080);
  // const Int_t Nx = 1;
  // const Int_t Ny = 3;
  // CanvasPartition(c, Nx, Ny);
  // TPad *pad[Nx][Ny];

  // for (Int_t i=0;i<Nx;i++) {
  //   for (Int_t j=0;j<Ny;j++) {
  //     c->cd(0);
  //     // Get the pads previously created.
  //     char pname[16];
  //     sprintf(pname,"pad_%i_%i",i,j);
  //     pad[i][j] = (TPad*) gROOT->FindObject(pname);
  //     pad[i][j]->Draw();
  //     pad[i][j]->SetFillStyle(4000);
  //     pad[i][j]->SetFrameFillStyle(4000);
  //     pad[i][j]->cd();
  //   }
  // }

  // pad[0][0]->cd();
  auto palette1 = new TPaletteAxis(132,0,142.5,1,h[0]);
  h[0]->GetListOfFunctions()->Add(palette1);

  c->Divide(1,2,0,0.15);
  c->cd(1);
  h[0]->Draw("colz");
  c->cd(2);
  h[1]->Draw("colz");


}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         auto name = TString::Format("pad_%d_%d",i,j);
         auto pad = (TPad*) C->FindObject(name.Data());
         if (pad) delete pad;
         pad = new TPad(name.Data(),"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}

