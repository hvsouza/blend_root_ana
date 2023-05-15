/*
======================================================================================
 *
 *       Filename:  calc.cpp
 *
 *    Description: Quick code for testin cmake stkills 
 *
 *        Version:  1.0
 *        Created:  05/15/2023 02:45:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */


#include <iostream>
#include <iostream>
#include "./calclib.cpp"
/* #include "./Event.cpp" */

#include "TGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRootCanvas.h"


int main(int argc, char **argv){
    std::cout << "Some code for testing" << std::endl;

    std::cout << division(5,2) << std::endl;
    
    std::cout << "plotting" << std::endl;
    /* plot_some_graph(argc, argv); */

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
    std::cout << "plotted" << std::endl;
    return 0;
}


