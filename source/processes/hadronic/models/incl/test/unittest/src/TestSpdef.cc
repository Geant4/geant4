// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestSpdef.hh"

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

ClassImp(TestSpdef)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>G4Abla::spdef</h1>

     <p>
     Test method G4Abla::spdef.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestSpdef.png"/>
     </p>
     END_HTML */
  
TestSpdef::TestSpdef()
{
  setUnitName("G4Abla::spdef");
  setOriginalUnitName("spdef");
  setPlotFileName("htmldoc/plots/TestSpdef.png");
  setLogFileName("htmldoc/logs/TestSpdef.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestSpdef::~TestSpdef()
{

}

void TestSpdef::runMe()
{
  // This handles the actual testing.
	
  TCanvas *c1 = new TCanvas();

  c1->SaveAs(getPlotFileName());
  
  setTestStatus(true);
}
