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
// $Id: TestFm2.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestFm2.hh"

#include "InclAblaTestDifferencePlotter.hh"

#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"

//#include "InclAblaFunctions.hh"
#include "functionwrapper.hh"
#include "commonwrapper.hh"

#ifndef __CINT__
#include "G4Incl.hh"
#include "G4InclDataDefs.hh"
#endif

#include <iostream>

ClassImp(TestFm2)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test fm2</h1>

     <p>
     This class tests method G4Incl::fm2
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestFm2.png"/>
     </p>
     END_HTML */
  
TestFm2::TestFm2()
{
  setUnitName("G4Incl::fm2");
  setOriginalUnitName("fm2");
  setPlotFileName("htmldoc/plots/TestFm2.png");
  setLogFileName("htmldoc/logs/TestFm2.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestFm2::~TestFm2()
{

}

void TestFm2::runMe()
{
  // This handles the actual testing.
	
  G4Incl *incl = new G4Incl();

  Float_t integral = 0.0;

  // Set the error marginal:
  const Float_t errorMarginal = 1e-3;

  const Int_t points = 20;
  Int_t j;
  Float_t r[points];
  Float_t cpp_f_r[points];
  Float_t fort_f_r[points];
  Float_t relativeDifference[points];

  for(Int_t i = 0; i < points; i = i + 1) {
    r[i] = i;
    cpp_f_r[i] = incl->fm2(i);
    fort_f_r[i] = fm2_(&i);
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    integral = integral + TMath::Abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, r, cpp_f_r);
  TGraph *fort_graph = new TGraph(points, r, fort_f_r);

  fort_graph->GetXaxis()->SetTitle("r");
  fort_graph->GetYaxis()->SetTitle("fm2(r)");
  fort_graph->SetTitle("Function fm2");
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

  // If the integral is smaller than our errorMarginal the test is passed.
  if(integral < errorMarginal) {
    // Passed the test
    setTestStatus(true);
  }
  else {
    // Failed the test
    setTestStatus(false);
  }
}
