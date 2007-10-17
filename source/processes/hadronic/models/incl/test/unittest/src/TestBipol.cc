// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestBipol.hh"

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

ClassImp(TestBipol)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Abla::bipol</h1>

     <p>
     Testing routine G4Abla::bipol.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestBipol.png"/>
     </p>
     <p>
     Log: <a href="logs/TestBipolCpp.txt">log for G4Abla::bipol</a> and for <a href="logs/TestBipolFortran.txt">FORTRAN bipol</a>
     </p>
     END_HTML */
  
TestBipol::TestBipol()
{
  setUnitName("G4Abla::bipol");
  setOriginalUnitName("bipol");
  setPlotFileName("htmldoc/plots/TestBipol.png");
  setLogFileName("htmldoc/logs/TestBipol.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestBipol::~TestBipol()
{

}

void TestBipol::runMe()
{
  // This handles the actual testing.
	
  double diff = 0.0;
  int rounds = 10;
  const double yStep = 0.1;
  double y = 0.0;
  int iflag = 0;
  double result, resultF;

  const double errorMarginal = 1e-6;
  
  G4Abla *abla = new G4Abla();

  ofstream log("htmldoc/logs/TestBipolCpp.txt");
  ofstream logF("htmldoc/logs/TestBipolFortran.txt");

  for(int i = 0; i < rounds; i++) {
    y = double(i)*yStep;
    result = abla->bipol(iflag, y);
    resultF = bipol_(&iflag, &y);
    log << "iflag = " << iflag << " y = " << y << " bipol = " << result << endl;
    logF << "iflag = " << iflag << " y = " << y << " bipol = " << resultF << endl;    
    diff = diff + fabs(result - resultF);
  }

  log.close();
  logF.close();
  
  TCanvas *c1 = new TCanvas();

  c1->SaveAs(getPlotFileName());

  if(diff < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
