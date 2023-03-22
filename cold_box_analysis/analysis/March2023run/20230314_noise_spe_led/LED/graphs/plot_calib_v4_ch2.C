
void plot_calib_v4_ch2(){

  // TFile *f1 = new TFile("./pe_vs_led_v4.root","READ");
  // TFile *f2 = new TFile("./pe_vs_led_v4_attenuator.root","READ");

  TFile *f1 = new TFile("./pe_vs_led_v4_Ch2.root","READ");
  TFile *f2 = new TFile("./pe_vs_led_v4_Ch2_attenuator.root","READ");

  TCanvas *c2 = (TCanvas*)f1->Get("c2");
  TGraphErrors *gn = (TGraphErrors*)c2->GetPrimitive("Graph")->Clone("gn");
  delete c2;
  f2->cd();
  c2 = (TCanvas*)f2->Get("c2");
  TGraphErrors *gatt = (TGraphErrors*)c2->GetPrimitive("Graph")->Clone("gatt");
  TCanvas *c = new TCanvas();
  gn->Draw("AP");
  gatt->SetMarkerColor(kRed);
  gatt->SetLineColor(kRed);
  gatt->SetMarkerSize(0.7);
  gn->SetMarkerSize(0.7);
  gatt->SetTitle("V4 Ch2 Attenuated (5 dB)");
  gn->SetTitle("V4 Ch2");
  gatt->Draw("SAME P");
  c->BuildLegend();
}
