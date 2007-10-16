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
// $Id: TestDerivWsax.cc,v 1.3 2007-10-16 20:37:44 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestDerivWsax.hh"

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

ClassImp(TestDerivWsax)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test Derivwsax</h1>

     <p>
     This class tests method G4Incl::derivWsax.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestDerivWsax.png"/>
     </p>
     END_HTML */
  
  TestDerivWsax::TestDerivWsax()
{
  setUnitName("G4incl::derivWsax");
  setOriginalUnitName("derivwsax");
  setPlotFileName("htmldoc/plots/TestDerivWsax.png");
  setLogFileName("htmldoc/logs/TestDerivWsax.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestDerivWsax::~TestDerivWsax()
{

}

void TestDerivWsax::runMe()
{
  // This handles the actual testing.

  // Create a new G4Incl instance:
  G4Incl *incl = new G4Incl();
  
  // Set starting values:
  Float_t r0 = 0.1;
  Float_t adif = 1.0;
  Float_t integral = 0.0;

  // Set the error marginal:
  const Float_t errorMarginal = 1e-6;

  // Set up arrays for data collection:
  const Int_t points = 50;
  const Float_t r_step = 0.3;
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

  for(Int_t i = 1; i < points; i++) {
    r[i] = r[i-1] + r_step;	
    cpp_f_r[i] = incl->derivWsax(r[i]);
    fort_f_r[i] = derivwsax_(&r[i]);
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    if(relativeDifference[i] > errorMarginal) {
      std::cout <<"C++: " << cpp_f_r[i] << "\t FORTRAN: " << fort_f_r[i]  << " Relative difference: " << relativeDifference[i] << std::endl;
    }
    integral = integral + abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, r, cpp_f_r);
  TGraph *fort_graph = new TGraph(points, r, fort_f_r);

  fort_graph->GetXaxis()->SetTitle("r");
  fort_graph->GetYaxis()->SetTitle("wsax(r)");
  fort_graph->SetTitle("Function derivwsax");
  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(points, r, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("r");
  diff_graph->GetYaxis()->SetTitle("Relative difference (%)");
  diff_graph->Draw("ap");

  c1->SaveAs(getPlotFileName());

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
