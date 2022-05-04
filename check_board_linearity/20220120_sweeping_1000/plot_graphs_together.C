
class DATA{

public:
  Double_t vnom;
  Double_t vin;
  Double_t vout;
  Double_t vinstd;
  Double_t voutstd;
  Double_t areain;
  Double_t areainstd;
  Double_t areaout;
  Double_t areaoutstd;
  Double_t offset;
  Double_t offsetstd;
  Double_t dataset;
  Double_t basein;
  Double_t baseinstd;
  Double_t baseout;
  Double_t baseoutstd;
};

void plot_graphs_together(){

  TFile *f = new TFile("cold_data.root","READ");
  TTree *t1 = (TTree*)f->Get("t1");
  Int_t n = t1->GetEntries();

  TFile *f2 = new TFile("warm_data.root","READ");
  TTree *t2 = (TTree*)f2->Get("t1");
  Int_t n2 = t2->GetEntries();


  DATA cold;
  DATA warm;
  // you dont need to read everybody
  t1->SetBranchAddress("vnom",&cold.vnom);
  t1->SetBranchAddress("vin",&cold.vin);
  t1->SetBranchAddress("vout",&cold.vout);
  t1->SetBranchAddress("vinstd",&cold.vinstd);
  t1->SetBranchAddress("voutstd",&cold.voutstd);
  t1->SetBranchAddress("areain",&cold.areain);
  t1->SetBranchAddress("areainstd",&cold.areainstd);
  t1->SetBranchAddress("areaout",&cold.areaout);
  t1->SetBranchAddress("areaoutstd",&cold.areaoutstd);
  t1->SetBranchAddress("offset",&cold.offset);
  t1->SetBranchAddress("offsetstd",&cold.offsetstd);
  t1->SetBranchAddress("dataset",&cold.dataset);
  t1->SetBranchAddress("basein",&cold.basein);
  t1->SetBranchAddress("baseinstd",&cold.baseinstd);
  t1->SetBranchAddress("baseout",&cold.baseout);
  t1->SetBranchAddress("baseoutstd",&cold.baseoutstd);

  t2->SetBranchAddress("vnom",&warm.vnom);
  t2->SetBranchAddress("vin",&warm.vin);
  t2->SetBranchAddress("vout",&warm.vout);
  t2->SetBranchAddress("vinstd",&warm.vinstd);
  t2->SetBranchAddress("voutstd",&warm.voutstd);
  t2->SetBranchAddress("areain",&warm.areain);
  t2->SetBranchAddress("areainstd",&warm.areainstd);
  t2->SetBranchAddress("areaout",&warm.areaout);
  t2->SetBranchAddress("areaoutstd",&warm.areaoutstd);
  t2->SetBranchAddress("offset",&warm.offset);
  t2->SetBranchAddress("offsetstd",&warm.offsetstd);
  t2->SetBranchAddress("dataset",&warm.dataset);
  t2->SetBranchAddress("basein",&warm.basein);
  t2->SetBranchAddress("baseinstd",&warm.baseinstd);
  t2->SetBranchAddress("baseout",&warm.baseout);
  t2->SetBranchAddress("baseoutstd",&warm.baseoutstd);

  
  Double_t vnom[n],vin[n],vout[n],vinstd[n],voutstd[n],areain[n],areainstd[n],areaout[n],areaoutstd[n],offset[n],offsetstd[n],dataset[n],basein[n],baseinstd[n],baseout[n],baseoutstd[n];
  Double_t vnom2[n2],vin2[n2],vout2[n2],vinstd2[n2],voutstd2[n2],areain2[n2],areainstd2[n2],areaout2[n2],areaoutstd2[n2],offset2[n2],offsetstd2[n2],dataset2[n2],basein2[n2],baseinstd2[n2],baseout2[n2],baseoutstd2[n2];

 
  
  Double_t corrected_areaout[n];
  Double_t corrected_areaout2[n2];
  Double_t corrected_areaoutstd[n];
  Double_t corrected_areaoutstd2[n2];

  Double_t const_out_baseline = 0;
  Double_t const_out_baselinestd = 0;

  Double_t convert_to_micro = 1000;
  Double_t charge_spe_cold = 7.78/convert_to_micro;
  Double_t charge_spe_warm = 7.8/convert_to_micro;

  for(Int_t i = 0; i<n; i++){
    t1->GetEntry(i);
    if(cold.vnom == 0){
      const_out_baseline = cold.baseout;
      const_out_baselinestd = cold.baseoutstd;
    }

    vnom[i] = cold.vnom;
    vin[i] = cold.vin*convert_to_micro;
    vout[i] = cold.vout;
    vinstd[i] = cold.vinstd*convert_to_micro;
    voutstd[i] = cold.voutstd;
    areain[i] = (cold.areain - cold.basein)/convert_to_micro/charge_spe_cold;
    areainstd[i] = (sqrt(pow(cold.areainstd,2)+pow(cold.baseinstd,2)))/convert_to_micro/charge_spe_cold;
    areaout[i] = cold.areaout-const_out_baseline;
    areaoutstd[i] = sqrt(pow(cold.areaoutstd,2)+pow(const_out_baselinestd,2)); ;
    corrected_areaout[i] = cold.areaout-cold.baseout;
//     cout <<  areaout[i] << " " << cold.areaout << " " << const_out_baseline << endl;
    corrected_areaoutstd[i] = sqrt(pow(cold.areaoutstd,2)+pow(cold.baseoutstd,2));
    offset[i] = cold.offset;
    offsetstd[i] = cold.offsetstd;
    dataset[i] = cold.dataset;

  }

  for(Int_t i = 0; i<n2; i++){
    t2->GetEntry(i);
    if(warm.vnom == 0){
      const_out_baseline = warm.baseout;
      const_out_baselinestd = warm.baseoutstd;
    }

    vnom2[i] = warm.vnom;
    vin2[i] = warm.vin*convert_to_micro;
    vout2[i] = warm.vout;
    vinstd2[i] = warm.vinstd*convert_to_micro;
    voutstd2[i] = warm.voutstd;
    areain2[i] = (warm.areain - warm.basein)/convert_to_micro/charge_spe_warm;
    areainstd2[i] = (sqrt(pow(warm.areainstd,2)+pow(warm.baseinstd,2)))/convert_to_micro/charge_spe_warm;
    areaout2[i] = warm.areaout-const_out_baseline;
    areaoutstd2[i] = sqrt(pow(warm.areaoutstd,2)+pow(const_out_baselinestd,2)); ;
    corrected_areaout2[i] = warm.areaout-warm.baseout;
//     cout <<  areaout2[i] << " " << warm.areaout << " " << const_out_baseline << endl;
    corrected_areaoutstd2[i] = sqrt(pow(warm.areaoutstd,2)+pow(warm.baseoutstd,2));
    offset2[i] = warm.offset;
    offsetstd2[i] = warm.offsetstd;
    dataset2[i] = warm.dataset;

  }




  // TCanvas *c1 = new TCanvas();
  // gStyle->SetOptFit(1);
  // TGraphErrors *g = new TGraphErrors(n,areain,areaout,areainstd,areaoutstd);
  // g->Draw("AP");
  // g->SetMarkerStyle(7);
  // g->SetNameTitle("area output vs input","area output vs input");
  // g->GetYaxis()->SetTitle("Area output (nV.s)");
  // g->GetXaxis()->SetTitle("Area input (nV.s)");

  // g->Fit("pol1");

  // TCanvas *c2 = new TCanvas();
  // TGraphErrors *g2 = new TGraphErrors(n,vin,vout,vinstd,voutstd);
  // g2->Draw("AP");
  // g2->SetMarkerStyle(7);
  // g2->SetNameTitle("amplitude output vs input","amplitude output vs input");
  // g2->GetYaxis()->SetTitle("Amplitude output (mV)");
  // g2->GetXaxis()->SetTitle("Amplitude input (#muV)");

  // g2->Fit("pol1");

  // TCanvas *c3 = new TCanvas();
  // TGraphErrors *g3 = new TGraphErrors(n,vnom,offset,0,offsetstd);
  // g3->Draw("AP");
  // g3->SetMarkerStyle(7);
  // g3->SetNameTitle("offset","offset");
  // g3->GetYaxis()->SetTitle("Offset amplitude (mV)");
  // g3->GetXaxis()->SetTitle("Nominal pulse amplitude (mV)");

  TCanvas *c4 = new TCanvas();
  TGraphErrors *g4 = new TGraphErrors(n,areain,corrected_areaout,areainstd,corrected_areaoutstd);
  TGraphErrors *g4warm = new TGraphErrors(n2,areain2,corrected_areaout2,areainstd2,corrected_areaoutstd2);
  g4->Draw("AP");
  g4warm->Draw("SAME P");
  g4->SetMarkerStyle(7);
  g4warm->SetMarkerStyle(7);
  g4->SetMarkerColor(kBlue);
  g4warm->SetMarkerColor(kRed);
  g4->SetLineColor(kBlue);
  g4warm->SetLineColor(kRed);
  g4->SetNameTitle("area output (corrected) vs input","area output (corrected) vs input");
  g4->GetYaxis()->SetTitle("Corrected area output (nV.s)");
  g4->GetXaxis()->SetTitle("Input photo-electrons");

  TF1 *flcold = new TF1("flcold","pol1",0,6);
  TF1 *flwarm = new TF1("flwarm","pol1",0,6);
  flcold->SetParameters(1,1);flcold->SetLineColor(kBlue);
  flwarm->SetParameters(1,1);flwarm->SetLineColor(kRed);
  g4->Fit("flcold");
  g4warm->Fit("flwarm");





}
