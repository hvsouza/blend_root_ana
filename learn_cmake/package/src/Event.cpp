/*
 * =====================================================================================
 *
 *       Filename:  Event.cpp
 *
 *    Description:  Some root plotting 
 *
 *        Version:  1.0
 *        Created:  05/15/2023 04:51:42 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */



#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include "TApplication.h"
#include "TRootCanvas.h"


void plot_some_graph(int argc, char **argv){

    TApplication app("app", &argc, argv);
    double x[5];
    double y[5];
    for(int i = 0; i < 5; i++){
        x[i] = i;
        y[i] = i*i;
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();

    TGraph *g = new TGraph(5, x, y);
    g->Draw("ALP*");
    c1->Modified(); c1->Update();

    TRootCanvas *rc = (TRootCanvas *)c1->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    app.Run();

    return;
}
