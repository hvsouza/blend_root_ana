// will just run example, so suggest if multiple files are needed
// use this to make all the rest ..
//root -b -q 'example.C()'
//help functions for making graphs and tex??
//root style file -> parent!
//watch out for directory structure
// #include "../rootstyle.C"


void example(const char dir[]=".")
{
  //make a plot, make some tex?
  // rootstyle();
  
  TCanvas *c1=new TCanvas("c1","my title",600,600);
  TGraph *g= new TGraph();
  g->SetPoint(0,10,10);
  g->SetPoint(1,20,20);
  g->SetPoint(2,30,30);

  g->Draw("Ap");
  g->GetXaxis()->SetTitle("xblah");
  g->GetYaxis()->SetTitle("yblah");

  c1->SaveAs(Form("%s/example_macro.png",dir)); // which format?
  std::ofstream fout(Form("%s/example_macro.tex",dir));
  fout<<"Some tex to go with example..."<<std::endl;
  //includegraphics ?
  fout.close();
  

}
#define A 0.8  
#define k 0.0347 //kV/MeV  stopping muons..
#define dEdx 2 //MIPs
double modBirks(double E)
{

  return (A/ (1 + (k/E)*dEdx));
}

Double_t fvel(Double_t *x, Double_t *par)
  {
    // A prototype liquid Argon Time Projection Chamber for the study of UV lasermulti-photonic ionization Article in Journal of Instrumentation · June 2009 DOI: 10.1088/1748-0221/4/07/P07011 · Source: arXi
    //vecoity in mm/microsecond
							    
  //Walkowiak measurements? not good below 0.5 kV/cm??
  double T=par[0];
  double E=x[0];
  double T0=90.371;//?
  double P1 = -0.01481; //+/- 0.00095 K^-1
  double P2 = -0.0075; //+/- 0.0028 K^-1
  double P3 = 0.141; //+/- 0.023 (kV/cm)^-1
  double P4 = 12.4; //+/- 2.7 kV/cm
  double P5 = 1.627; //+/- 0.078 (kV/cm)^P6
  double P6 = 0.317; //+/- 0.021
  
  return (P1*(T-T0) + 1)*(P3*E*TMath::Log(1+ (P4/E)) + P5*pow(E,P6)) + P2*(T-T0);

}
void lar_vel(const char dir[]=".")
{
  //make a plot, make some tex?
  // rootstyle();

  double Eliq,T;
  T=88;
  
  TCanvas *c1=new TCanvas("c1","my title",600,600);
  TGraph *g= new TGraph();
  Eliq=0;
  double Eliqstep=5.0/100;
  for(int i=0;i<100;i++)
    {
      Eliq = Eliq + Eliqstep;
      g->SetPoint(g->GetN(), Eliq,fvel(&Eliq, &T) );
    }

  g->Draw("AC");
  g->GetXaxis()->SetTitle("Electric Field Strength in Liquid Argon (kV/cm) ");
  g->GetYaxis()->SetTitle("Electron drift speed mm/#mus");

  c1->SaveAs(Form("%s/electrondriftspeed.png",dir)); // which format?
  //  std::ofstream fout(Form("%s/example_macro.tex",dir));
  // fout<<"Some tex to go with example..."<<std::endl;
  //includegraphics ?
  // fout.close();

 }
