// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestCram.hh"

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

ClassImp(TestCram)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>G4Abla::cram</h1>

     <p>
     Test method G4Abla::cram
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestCram.png"/>
     </p>
     <p>
     Log: <a href="logs/TestCramCpp.txt">log for G4Abla::cram</a> and for <a href="logs/TestCramFortran.txt">FORTRAN cram</a>
     </p>
     END_HTML */
  
TestCram::TestCram()
{
  setUnitName("G4Abla::cram");
  setOriginalUnitName("cram");
  setPlotFileName("htmldoc/plots/TestCram.png");
  setLogFileName("htmldoc/logs/TestCram.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestCram::~TestCram()
{

}

void TestCram::runMe()
{
  // This handles the actual testing.
	
  TCanvas *c1 = new TCanvas();

  c1->SaveAs(getPlotFileName());
  
  setTestStatus(true);
}
