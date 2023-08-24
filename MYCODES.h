// ________________________________________ //
// Author: Henrique Souza
// Filename: MYCODES.h
// Created: 2021
// ________________________________________ //
#ifndef memorydepth
#define memorydepth 5000
#endif

#ifndef memorydepth_sample
#define memorydepth_sample memorydepth
#endif

#ifndef MYCODES
#define MYCODES

#include <fstream>
#include <iostream>
#include "Riostream.h"
#include <time.h>       // time_t, struct tm, difftime, time, mktime
#include <cmath>
#include <numeric>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <cstring> // include the <cstring> header file for memset()


#include "TROOT.h"
#include "TLatex.h"
#include <TMinuit.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TLine.h>
#include <TTimeStamp.h>
#include <TComplex.h>
#include <TVirtualFFT.h>
#include <TMatrixD.h>
#include "TSystem.h"
#include <TEventList.h>
#include <TPolyMarker.h>
#include <TSpectrum.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TMarker.h>
#include <TKey.h>


// #include "DUNEStyle.h"
// // use with:
// //
// // #define DUNESTYLE_ENABLE_AUTOMATICALLY 0
// // dunestyle::SetDuneStyle();

using namespace std;

#include "ADC_DATA.h"
#include "old_class.C"
#include "timeReader.C"
#include "readingCodes.C"
#include "wiener_filter.C"
#include "analyzer.C"
#include "analyzer_plots.C"
#include "calibrationCodes.C"
#include "old_search_spe.C"
#include "alpha_analysis.C"
#include "samples.C"

// #include "calibrationCodes_CT_validation.C"


#endif


