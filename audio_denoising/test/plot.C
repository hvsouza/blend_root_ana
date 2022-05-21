void plot(){

  TCanvas *cv = new TCanvas("cv","cv");
  TGraph *gv = new TGraph("v.txt");
  gv->Draw("ALP");

  TCanvas *cu = new TCanvas("cu","cu");
  TGraph *gu = new TGraph("u.txt");
  gu->Draw("ALP");

  TCanvas *cs = new TCanvas("cs","cs");
  TGraph *gs = new TGraph("s.txt");
  gs->Draw("ALP");

  TCanvas *cx = new TCanvas("cx","cx");
  TGraph *gx = new TGraph("x.txt");
  gx->Draw("ALP");

  TCanvas *cy = new TCanvas("cy","cy");
  TGraph *gy = new TGraph("y.txt");
  gy->Draw("ALP");

}
