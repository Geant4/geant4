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
// $Id: TestRadius.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestRadius.hh"

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

#include <iostream>

ClassImp(TestRadius)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test radius</h1>

     <p>
     This class tests method G4Incl::radius.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestRadius.png"/>
     </p>
     END_HTML */
  
TestRadius::TestRadius()
{
  setUnitName("G4Incl::radius");
  setOriginalUnitName("radi_us");
  setPlotFileName("htmldoc/plots/TestRadius.png");
  setLogFileName("htmldoc/logs/TestRadius.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestRadius::~TestRadius()
{

}

void TestRadius::runMe()
{
  // This handles the actual testing.

  const double errorMarginal = 1e-3;
  
  G4Incl *incl = new G4Incl();

  double minA = 28.0;
  const int nuclei = 50;

  double A = minA;

  double Ai[nuclei];
  double cppRadius[nuclei];
  double fortRadius[nuclei];

  double relativeDifference[nuclei];
  double integral = 0.0;
  
  for(Int_t i = 0; i < nuclei; i++) {
    Ai[i] = A;
    cppRadius[i] = incl->radius(Ai[i]);
    fortRadius[i] = radi_us__(&Ai[i]);
    //    cout << cppRadius[i] << " \t " << fortCX << endl;
    relativeDifference[i] = 100.0*(cppRadius[i] - fortRadius[i])/fortRadius[i];
    integral = integral + abs(cppRadius[i] - fortRadius[i]);
    A = A + 1.0;
  }
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(nuclei, Ai, cppRadius);
  TGraph *fort_graph = new TGraph(nuclei, Ai, fortRadius);

  fort_graph->GetXaxis()->SetTitle("A");
  fort_graph->GetYaxis()->SetTitle("sig_reac(bullet=proton, E=50MeV, A)");
  fort_graph->SetTitle("Function sig_reac");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(nuclei, Ai, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("r");
  diff_graph->GetYaxis()->SetTitle("Relative difference (%)");

  diff_graph->Draw("ap");

  c1->SaveAs(getPlotFileName());

  // Clean up...
  delete c1;
  delete cpp_graph;
  delete fort_graph;
  delete diff_graph;

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
