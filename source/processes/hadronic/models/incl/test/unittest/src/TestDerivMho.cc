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
// $Id: TestDerivMho.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestDerivMho.hh"

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

ClassImp(TestDerivMho)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test derivmho</h1>

     <p>
     This class tests method G4Incl::derivMho.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestDerivMho.png"/>
     </p>
     END_HTML */
  
  TestDerivMho::TestDerivMho()
{
  setUnitName("G4Incl::derivMho");
  setOriginalUnitName("derivmho");
  setPlotFileName("htmldoc/plots/TestDerivMho.png");
  setLogFileName("htmldoc/logs/TestDerivMho.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestDerivMho::~TestDerivMho()
{

}

void TestDerivMho::runMe()
{
  // This handles the actual testing.

  // New G4Incl instance
  G4Incl *incl = new G4Incl();

  // Error marginal
  const Float_t errorMarginal = 1e-9;
  
  // Set starting values:
  Float_t r0 = 0.1;
  Float_t adif = 1.0;
  Float_t integral = 0.0;

  const Int_t points = 50;
  const Float_t r_step = 0.2;
  Float_t r[points];
  r[0] = 0.0;
  Float_t cpp_f_r[points];
  Float_t fort_f_r[points];
  Float_t relativeDifference[points];

  // C++:
  G4Dton *dton = (G4Dton*) malloc(sizeof(G4Dton));
  G4Ws *ws = (G4Ws*) malloc(sizeof(G4Ws));
  ws->r0 = r0;
  ws->adif = adif;

  incl->setDtonData(dton);
  incl->setWsData(ws);

  // FORTRAN:
  gWs->r0 = r0;
  gWs->adif = adif;

  // Call C++ and FORTRAN functions:
  for(Int_t i = 1; i < points; i++) {
    r[i] = r[i-1] + r_step;	
    cpp_f_r[i] = incl->derivMho(r[i]);
    fort_f_r[i] = derivmho_(&r[i]);
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    integral = integral + abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, r, cpp_f_r);
  TGraph *fort_graph = new TGraph(points, r, fort_f_r);

  fort_graph->GetXaxis()->SetTitle("r");
  fort_graph->GetYaxis()->SetTitle("derivmho(r)");
  fort_graph->SetTitle("Function derivmho");
  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(points, r, relativeDifference);
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

  // The integral over the difference of the bins of the histograms should be 0
  // if both C++ and FORTRAN results are the same.

  if(integral < errorMarginal) {
    // Passed the test
    setTestStatus(true);
  }
  else {
    // Failed the test
    setTestStatus(false);
  }
}
