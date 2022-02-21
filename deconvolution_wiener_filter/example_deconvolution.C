/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : example_deconvolution
 * @created     : mercoledì gen 12, 2022 13:31:21 CET
 */

#include <stdio.h>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TComplex.h"

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
  T h = (b - a) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}

// compute the integral of a TGraph object
double g_integral(TGraph* g, double x0, double x1) {
  auto fc_gintegral = [g](double *x, double* p) {
    return p[0]*g->Eval(x[0]);
  };
  TF1 f("f", fc_gintegral, x0, x1, 1);
  f.SetNpx(g->GetN());
  f.SetParameter(0, 1.);
  return f.Integral(x0, x1,1e-4);
}

// scale all the values of a TGraph
void g_scale(TGraph* g, double c) {
  for (int i=0; i<g->GetN(); i++) {
    g->GetY()[i] *= c;
  }
  return;
}

TGraph* build_noise_spectral_density(double rms, int nwindow, int nsample, double t1, double t0) {
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;

  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;

  int nsample_fft = 0.5*nsample+1;
  TGraph* g_noise_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++) 
    g_noise_spectral_density->SetPoint(j, j*0.25*TMath::InvPi()/(t1-t0), 0.);
  double ymin = 1e-7; 
  double ymax = 1e-3;
  int    nbinsy = 100;
  
  std::vector<double> ybin_exponent = 
    linspace(TMath::Log10(ymin), TMath::Log10(ymax), nbinsy);
  std::vector<double> ybins(100, 0);
  for (int iy=0; iy<100; iy++) {
    ybins[iy] = TMath::Power(10., ybin_exponent[iy]);
  }

  TH2D* h2_spectral_density = new TH2D("h2_spectral_density", 
      Form("%s;%s;%s", "Noise spectral density", "Frequency [A.U.]", "Amplitude [A.U.]"),
      nsample_fft, 
      g_noise_spectral_density->GetXaxis()->GetXmin(),
      g_noise_spectral_density->GetXaxis()->GetXmax(),
      nbinsy-1, &ybins.at(0));
  h2_spectral_density->GetXaxis()->CenterTitle();
  h2_spectral_density->GetYaxis()->CenterTitle();

  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");
  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) {
      xn[ip] = gRandom->Gaus(0, rms);
    }

    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j]) * c_scale;
      h2_spectral_density->Fill(j*0.25*TMath::InvPi()/(t1-t0), xN[j].Rho2());
      g_noise_spectral_density->GetY()[j] += (xN[j].Rho2()*scale);
    }
  }

  TCanvas* cNoise = new TCanvas("cNoise", "Noise spectral density", 100, 100, 800, 600);
  cNoise->SetLogy(1);
  cNoise->SetTicks(1, 1);
  cNoise->SetGrid(1, 1);
  h2_spectral_density->Draw("col");
  g_noise_spectral_density->SetLineColor(kGray+2);
  g_noise_spectral_density->SetLineWidth(2);
  g_noise_spectral_density->Draw("l");
  return g_noise_spectral_density;
}

int example_deconvolution()
{
  gStyle->SetPalette(kSunset);
  //- - - - - - - - - - - - - - - - - - - - - - Retrieve spe template waveform

  TFile* template_file = new TFile("./spe_template_62_5MHz.root");
  TGraph* spe_template = (TGraph*)template_file->Get("spe_template_62_5MHz");

  //- - - - - - - - - - - - - - - - - - - - - - - - -Build a syntetic waveform

  double dt = 1.6e-02; // in us
  double t0 = 0.;      // in us
  double t1 = 16;      // in us
  const int nsample = 1000;
  const double tmpl_prepulse_us   = 1.1; // template pre-pulse time
  const int tmpl_prepulse_tick =  68; // template pre-pulse ticks:w

  gRandom->SetSeed(0);
  const int nph_true = gRandom->Poisson(200);   // true number of p.e. 
  std::vector<double> tph_true(nph_true, 0.); // true time of p.e. 
  std::vector<TLine*> lph_true(nph_true, 0);  // lines for p.e. display

  // Sample true p.e. time from and exp with τ=1.6 μs and shift 5 μs
  for (int i=0; i<nph_true; i++) {
    tph_true[i] = gRandom->Exp(2.6)+5;    //1.6
    lph_true[i] = new TLine(tph_true[i], 0., tph_true[i], 4.);
    lph_true[i]->SetLineColor(kBlack);
    lph_true[i]->SetLineWidth(2);
  }

  // prepare auxiliary arrays
  double t = t0;
  double xt[nsample] = {0}; // time array
  double xv[nsample] = {1}; // waveform array
  double xn[nsample] = {0}; // noise array
  double xh[nsample] = {0}; // impulse response function array (spe template)
  double xs[nsample] = {0}; // original signal (δ-like function)
  double* xy;               // deconvoluted signal

  // fill impulse response function array
  for (int ip=0; ip<spe_template->GetN()-tmpl_prepulse_tick; ip++) {
    xh[ip] = spe_template->GetY()[ip+68];
  }

  // build noise spectral density
  TGraph* gNoise_spectral_density = 
    build_noise_spectral_density(0.25, 100, nsample, t1, t0);

  // build waveform and fill the above arrays
  for (int i=0; i<nsample; i++) { 
    for (int j=0; j<nph_true; j++) {
      double tt = t+tmpl_prepulse_us -tph_true[j];
      if (tt>1 && tt<9) { // avoid evaluating the template outside its domain 
        double v = spe_template->Eval(tt, 0, "S");
        xv[i] += v;
      }
    }
    xs[i] = TMath::Gaus(t, 1.1, 0.01, true)/dt; // true signal shape
    xt[i] = t;
    xn[i] = gRandom->Gaus(0, 0.25); // white noise
    xv[i] += xn[i]; // waveform = signal + noise
    t+=dt;
  }

  // display waveform and noise
  TGraph* gv = new TGraph(nsample, xt, xv);
  TGraph* gn = new TGraph(nsample, xt, xn);
  TGraph* gs = new TGraph(nsample, xt, xs);

  gv->SetLineColor(kRed+1);
  gn->SetLineColor(kGray+1);
  gs->SetLineColor(kBlue+1);

  gStyle->SetOptTitle(0);
  TCanvas* cTime = new TCanvas();
  cTime->Divide(1, 2);
  cTime->cd(1);
  gv->Draw("awlx+");
  gv->SetNameTitle("gv", "Syntetic waveform");
  gv->GetXaxis()->SetTitle("Time [#mus]");
  gv->GetYaxis()->SetTitle("Amplitude [mV]");
  gv->GetXaxis()->CenterTitle();
  gv->GetYaxis()->CenterTitle();
  gv->GetYaxis()->SetTitleSize(0.06);
  gv->GetYaxis()->SetLabelSize(0.06);
  gv->GetXaxis()->SetTitleSize(0.06);
  gv->GetXaxis()->SetLabelSize(0.06);
  //gn->Draw("l");
  //gs->Draw("l");
  gPad->SetTopMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.01);
  gPad->SetTicks(1, 1); gPad->SetGrid(1, 1);

  for (const auto &line : lph_true) line->Draw("same");
  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - Perform FFT

  int nsample_ = nsample;
  // xV: FFT of waveform
  TComplex xV[nsample]; double xV_re[nsample]; double xV_im[nsample];
  // xH: FFT of spe response
  TComplex xH[nsample]; double xH_re[nsample]; double xH_im[nsample];
  // xS: FFT of original signal
  TComplex xS[nsample]; double xS_re[nsample]; double xS_im[nsample];
  // xY: FFT of the filtered signal
  TComplex xY[nsample]; double xY_re[nsample]; double xY_im[nsample];

  // Instance the FFT engine
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample_, "M R2C");
  // wavaform FFT
  fft->SetPoints(xv);
  fft->Transform();
  fft->GetPointsComplex(xV_re, xV_im);

  // spe template FFT
  fft->SetPoints(xh);
  fft->Transform();
  fft->GetPointsComplex(xH_re, xH_im);

  // Original signal FFT
  fft->SetPoints(xs);
  fft->Transform();
  fft->GetPointsComplex(xS_re, xS_im);

  // Fill FFT arrays and perform Wiener deconvolution 
  double c_scale = 1./nsample;
  double H2[nsample] = {1}; // spe response spectral density
  double N2[nsample] = {0}; // noise spectral density
  double S2[nsample] = {0}; // original signal spectral density
  double xf[nsample] = {0}; // frequency array
  TComplex G[nsample];
  for (int i=0; i<nsample*0.5+1; i++) {
    // fill FFT arrays
    xH[i] = TComplex(xH_re[i], xH_im[i]) * c_scale;
    xV[i] = TComplex(xV_re[i], xV_im[i]) * c_scale;
    xS[i] = TComplex(xS_re[i], xS_im[i]) * c_scale;
    // Compute spectral sensity
    H2[i] = xH[i].Rho2();
    N2[i] = gNoise_spectral_density->GetY()[i];
    S2[i] = xS[i].Rho2();
    // Compute Wiener filter
    G[i]  = TComplex::Conjugate(xH[i])*S2[i] / (H2[i]*S2[i] + N2[i]);
    xf[i] = i*0.25*TMath::InvPi()/(t1-t0);

    // Compute filtered signal
    xY[i] = G[i]*xV[i];
    xY_re[i] = xY[i].Re(); xY_im[i] = xY[i].Im();
  }
    
  // Display (normalized) spectral densities
  TGraph* gN2 = new TGraph(nsample*0.5+1, xf, N2);
  TGraph* gH2 = new TGraph(nsample*0.5+1, xf, H2);
  TGraph* gS2 = new TGraph(nsample*0.5+1, xf, S2);

  double gN2_int = g_integral(gN2, 0., xf[(int)(nsample*0.5)]);
  double gH2_int = g_integral(gH2, 0., xf[(int)(nsample*0.5)]);
  double gS2_int = g_integral(gS2, 0., xf[(int)(nsample*0.5)]);

  printf("N2 integral = %g\n", gN2_int);
  printf("H2 integral = %g\n", gH2_int);
  printf("S2 integral = %g\n", gS2_int);

  g_scale(gN2, 1./gN2_int);
  g_scale(gH2, 1./gH2_int);
  g_scale(gS2, 1./gS2_int);

  TCanvas* cPower = new TCanvas(); 
  cPower->SetLogy(1);
  cPower->cd();
  gN2->SetLineColor(kGray+1);
  gH2->SetLineColor(kRed+1);
  gS2->SetLineColor(kBlue+1);
  gH2->Draw("awl");
  gN2->Draw("l");
  gS2->Draw("l");

  //- - - - - - - - - - - - - - - - -  Backward transform of the filtered signal
  fft = TVirtualFFT::FFT(1, &nsample_, "M C2R");
  fft->SetPointsComplex(xY_re, xY_im);
  fft->Transform();
  xy = fft->GetPointsReal();
  for (int i=0; i<nsample; i++) xy[i] *= nsample;

  cTime->cd(2);
  TGraph* gy = new TGraph(nsample, xt, xy);
  gy->GetXaxis()->SetTitle("Time [#mus]");
  gy->GetYaxis()->SetTitle("Amplitude [A.U.]");
  gy->SetTitle("Deconvolved waveform");
  gy->GetXaxis()->CenterTitle();
  gy->GetYaxis()->CenterTitle();
  gy->GetYaxis()->SetTitleSize(0.06);
  gy->GetYaxis()->SetLabelSize(0.06);
  gy->GetXaxis()->SetTitleSize(0.06);
  gy->GetXaxis()->SetLabelSize(0.06);
  gy->SetLineColor(kBlue+1);
  gy->SetLineWidth(2);
  gy->Draw("awl");
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTicks(1, 1); gPad->SetGrid(1, 1);
  gPad->BuildLegend(0.5, 0.88, 0.88, 0.8, "", "l");

  return 0;
}

