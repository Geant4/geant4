// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestLpoly.hh"

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

ClassImp(TestLpoly)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>G4Abla::Lpoly</h1>

     <p>
     Legendre polynomial test.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestLpoly.png"/>
     </p>
     <p>
     Log: <a href="logs/TestLpolyCpp.txt">log for G4Abla::lpoly</a> and for <a href="logs/TestLpolyFortran.txt">FORTRAN lpoly</a>
     </p>
     END_HTML */
  
TestLpoly::TestLpoly()
{
  setUnitName("G4Abla::Lpoly");
  setOriginalUnitName("lpoly");
  setPlotFileName("htmldoc/plots/TestLpoly.png");
  setLogFileName("htmldoc/logs/TestLpoly.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestLpoly::~TestLpoly()
{

}

void TestLpoly::runMe()
{
  // This handles the actual testing.

  const double errorMarginal = 1e-9;
  
  const int terms = 9;
  int fterms = terms;
  double x, pl[terms], plfortran[terms];

  double xStep = 0.1;
  
  ofstream logOut("htmldoc/logs/TestLpolyCpp.txt");
  ofstream logOutF("htmldoc/logs/TestLpolyFortran.txt");  

  TCanvas *c1 = new TCanvas();

  G4Abla *abla = new G4Abla();
  G4double diff = 0.0;
  
  logOut <<"Legendre pol. terms: " << endl;
  for(int i = 0; i < 100; i++) {
    x = xStep*double(i);
    abla->lpoly(x, terms, pl);
    lpoly_(&x, &fterms, plfortran);
    for(int j = 0; j < terms; j++) {
      logOut << pl[j] << " \t ";
      logOutF << plfortran[j] << " \t ";
      diff = diff + fabs(pl[j] - plfortran[j]);
    }
    logOut << endl;
    logOutF << endl;
  }
  
  c1->SaveAs(getPlotFileName());

  logOut.close();
  logOutF.close();

  if(diff < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
