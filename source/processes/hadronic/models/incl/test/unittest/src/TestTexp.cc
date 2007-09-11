//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TestTexp.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestTexp.hh"

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

ClassImp(TestTexp)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test texp</h1>

     <p>
     This class tests method G4Incl::safeExp.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestTexp.png"/>
     </p>
     END_HTML */
  
TestTexp::TestTexp()
{
  setUnitName("G4Incl::safeExp");
  setOriginalUnitName("texp");
  setPlotFileName("htmldoc/plots/TestTexp.png");
  setLogFileName("htmldoc/logs/TestTexp.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestTexp::~TestTexp()
{

}

void TestTexp::runMe()
{
  // This handles the actual testing.
	
  G4Incl *incl = new G4Incl();

  const double errorMarginal = 1e-9;
  const int points = 10;
  double r[points];
  double cpp_f_r[points];
  double fort_f_r[points];
  double relativeDifference[points];
  double r_step = 1.0;
  double r_start = 0.0;
  double integral = 0.0;

  for(int i = 0; i < points; i++) {
    if(i == 0) {
      r[i] = r_start;
    }
    else {
      r[i] = r[i-1] + r_step;
    }

    cpp_f_r[i] = incl->safeExp(r[i]);
    fort_f_r[i] = texp_(&r[i]);
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    integral = integral + abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();

  c1->Divide(2,2);

  TGraph *cppgraph = new TGraph(points, r, cpp_f_r);
  TGraph *fortgraph = new TGraph(points, r, fort_f_r);
  TGraph *differencegraph = new TGraph(points, r, relativeDifference);

  fortgraph->GetXaxis()->SetTitle("x");
  fortgraph->GetYaxis()->SetTitle("flin2(x)");
  fortgraph->SetTitle("Function flin2");

  c1->cd(1);
  fortgraph->Draw("al");
  cppgraph->Draw("p,same");

  c1->cd(4);
  differencegraph->Draw("ap");
  differencegraph->SetTitle("Relative difference C++/FORTRAN");
  differencegraph->GetXaxis()->SetTitle("r");
  differencegraph->GetYaxis()->SetTitle("Relative difference (%)");

  c1->SaveAs(getPlotFileName());
  
  // Clean up...
  delete c1;
  delete cppgraph;
  delete fortgraph;
  delete differencegraph;

  // The integral over the difference of the bins of the histograms should be
  // smaller than errorMarginal if both C++ and FORTRAN results are the same.
  if(integral < errorMarginal) {
    // Passed the test
    setTestStatus(true);
  }
  else {
    // Failed the test
    setTestStatus(false);
  }
}
