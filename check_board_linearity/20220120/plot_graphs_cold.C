void plot_graphs_cold(){

  TFile *f = new TFile("cold_data.root","READ");
  TTree *t1 = (TTree*)f->Get("t1");

  Int_t n = t1->GetEntries();

  Double_t vnom[n],vin[n],vout[n],vinstd[n],voutstd[n],areain[n],areainstd[n],areaout[n],areaoutstd[n],offset[n],offsetstd[n],dataset[n],basein[n],baseinstd[n],baseout[n],baseoutstd[n];
  Double_t r_vnom,r_vin,r_vout,r_vinstd,r_voutstd,r_areain,r_areainstd,r_areaout,r_areaoutstd,r_offset,r_offsetstd,r_dataset,r_basein, r_baseinstd,r_baseout,r_baseoutstd;

  // you dont need to read everybody
  t1->SetBranchAddress("vnom",&r_vnom);
  t1->SetBranchAddress("vin",&r_vin);
  t1->SetBranchAddress("vout",&r_vout);
  t1->SetBranchAddress("vinstd",&r_vinstd);
  t1->SetBranchAddress("voutstd",&r_voutstd);
  t1->SetBranchAddress("areain",&r_areain);
  t1->SetBranchAddress("areainstd",&r_areainstd);
  t1->SetBranchAddress("areaout",&r_areaout);
  t1->SetBranchAddress("areaoutstd",&r_areaoutstd);
  t1->SetBranchAddress("offset",&r_offset);
  t1->SetBranchAddress("offsetstd",&r_offsetstd);
  t1->SetBranchAddress("dataset",&r_dataset);
  t1->SetBranchAddress("basein",&r_basein);
  t1->SetBranchAddress("baseinstd",&r_baseinstd);
  t1->SetBranchAddress("baseout",&r_baseout);
  t1->SetBranchAddress("baseoutstd",&r_baseoutstd);

  Double_t corrected_areaout[n];
  Double_t corrected_areaoutstd[n];

  Double_t const_out_baseline = 0;
  Double_t const_out_baselinestd = 0;

  Double_t convert_to_micro = 1000;

  for(Int_t i = 0; i<n; i++){
    t1->GetEntry(i);
    if(r_vnom == 0){
      const_out_baseline = r_baseout;
      const_out_baselinestd = r_baseoutstd;
    }

    vnom[i] = r_vnom;
    vin[i] = r_vin*convert_to_micro;
    vout[i] = r_vout;
    vinstd[i] = r_vinstd*convert_to_micro;
    voutstd[i] = r_voutstd;
    areain[i] = (r_areain - r_basein)/convert_to_micro;
    areainstd[i] = (sqrt(pow(r_areainstd,2)+pow(r_baseinstd,2)))/convert_to_micro;
    areaout[i] = r_areaout-const_out_baseline;
    areaoutstd[i] = sqrt(pow(r_areaoutstd,2)+pow(const_out_baselinestd,2)); ;
    corrected_areaout[i] = r_areaout-r_baseout;
//     cout <<  areaout[i] << " " << r_areaout << " " << const_out_baseline << endl;
    corrected_areaoutstd[i] = sqrt(pow(r_areaoutstd,2)+pow(r_baseoutstd,2));
    offset[i] = r_offset;
    offsetstd[i] = r_offsetstd;
    dataset[i] = r_dataset;

  }



  TCanvas *c1 = new TCanvas();
  gStyle->SetOptFit(1);
  TGraphErrors *g = new TGraphErrors(n,areain,areaout,areainstd,areaoutstd);
  g->Draw("AP");
  g->SetMarkerStyle(7);
  g->SetNameTitle("area output vs input","area output vs input");
  g->GetYaxis()->SetTitle("Area output (nV.s)");
  g->GetXaxis()->SetTitle("Area input (nV.s)");

  g->Fit("pol1");

  TCanvas *c2 = new TCanvas();
  TGraphErrors *g2 = new TGraphErrors(n,vin,vout,vinstd,voutstd);
  g2->Draw("AP");
  g2->SetMarkerStyle(7);
  g2->SetNameTitle("amplitude output vs input","amplitude output vs input");
  g2->GetYaxis()->SetTitle("Amplitude output (mV)");
  g2->GetXaxis()->SetTitle("Amplitude input (#muV)");

  g2->Fit("pol1");

  TCanvas *c3 = new TCanvas();
  TGraphErrors *g3 = new TGraphErrors(n,vnom,offset,0,offsetstd);
  g3->Draw("AP");
  g3->SetMarkerStyle(7);
  g3->SetNameTitle("offset","offset");
  g3->GetYaxis()->SetTitle("Offset amplitude (mV)");
  g3->GetXaxis()->SetTitle("Nominal pulse amplitude (mV)");

  TCanvas *c4 = new TCanvas();
  TGraphErrors *g4 = new TGraphErrors(n,areain,corrected_areaout,areainstd,corrected_areaoutstd);
  g4->Draw("AP");
  g4->SetMarkerStyle(7);
  g4->SetNameTitle("area output (corrected) vs input","area output (corrected) vs input");
  g4->GetYaxis()->SetTitle("Corrected area output (nV.s)");
  g4->GetXaxis()->SetTitle("Area input (nV.s)");

  g4->Fit("pol1");





}
