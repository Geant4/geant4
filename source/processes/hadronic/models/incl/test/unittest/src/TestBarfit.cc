// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestBarfit.hh"

#include "InclAblaTestDifferencePlotter.hh"

#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"

//#include "InclAblaFunctions.hh"
#include "functionwrapper.hh"
#include "commonwrapper.hh"

#ifndef __CINT__
#include "G4Incl.hh"
#include "G4InclDataDefs.hh"
#endif

#include "Riostream.h"

ClassImp(TestBarfit)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>G4Abla::barfit</h1>

     <p>
     Testing method G4Abla::barfit.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestBarfit.png"/>
     </p>
     <p>
     Log: <a href="logs/TestBarfitCpp.txt">log for G4Abla::lpoly</a> and for <a href="logs/TestBarfitFortran.txt">FORTRAN lpoly</a>
     </p>
     END_HTML */
  
TestBarfit::TestBarfit()
{
  setUnitName("G4Abla::barfit");
  setOriginalUnitName("barfit");
  setPlotFileName("htmldoc/plots/TestBarfit.png");
  setLogFileName("htmldoc/logs/TestBarfit.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestBarfit::~TestBarfit()
{

}

void TestBarfit::runMe()
{
  // This handles the actual testing.
  int ia = 150;
  int iz = 82;
  int il = 0;
  double sbfis = 0.0, segs = 0.0, selmax = 0.0;
  float sbfisF = 0.0, segsF = 0.0, selmaxF = 0.0;

  const int points = 100;
  double neutronnumber[points];
  double cppresult[points];
  double fortresult[points];
  
  G4Abla *abla = new G4Abla();

  ofstream log("htmldoc/logs/TestBarfitCpp.txt");
  ofstream logF("htmldoc/logs/TestBarfitFortran.txt");

  const double errorMarginal = 1e-3;
  double diff = 0.0;

  for(int i = 0; i < points; i++) {
    ia = ia + i;
    barfit_(&iz, &ia, &il, &sbfisF, &segsF, &selmaxF);
    abla->barfit(iz, ia, il, &sbfis, &segs, &selmax);
    log <<"ia = " << ia << " iz = " << iz << " sbfis = " << sbfis << " segs = " << segs << " selmax = " << selmax << endl;
    logF <<"ia = " << ia << " iz = " << iz << " sbfis = " << sbfisF << " segs = " << segsF << " selmax = " << selmaxF << endl; 
    neutronnumber[i] = ia-iz;
    cppresult[i] = sbfis;
    fortresult[i] = sbfisF;
    diff = diff + fabs(sbfis - sbfisF) + fabs(segs + segsF) + fabs(selmax - selmaxF);
  }
  
  TCanvas *c1 = new TCanvas();
  
  TGraph *cpp_graph = new TGraph(points, neutronnumber, cppresult);
  TGraph *fort_graph = new TGraph(points, neutronnumber, fortresult);
  fort_graph->Draw("al");

  cout <<"TestBarfit: plotting" << endl;
  
  //  fort_graph->SetMarkerStyle(2);

//   cpp_graph->Draw("al");
//   fort_graph->Draw("p, same");
  
  cout <<"TestBarfit: plot complete." << endl;
  //  c1->SaveAs(getPlotFileName());

  if(diff < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
    cout <<"Test for barfit failed. diff = " << diff << endl; 
  }
}
