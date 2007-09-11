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
// $Id: TestCoulombTransm.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestCoulombTransm.hh"

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

ClassImp(TestCoulombTransm)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::coulombTransm</h1>

     <p>
     Test method G4Incl::coulombTransm.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestCoulombTransm.png"/>
     </p>
     END_HTML */
  
TestCoulombTransm::TestCoulombTransm()
{
  setUnitName("G4Incl::coulombTransm");
  setOriginalUnitName("coulomb_transm");
  setPlotFileName("htmldoc/plots/TestCoulombTransm.png");
  setLogFileName("htmldoc/logs/TestCoulombTransm.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestCoulombTransm::~TestCoulombTransm()
{

}

void TestCoulombTransm::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();

  const double errorMarginal = 1e-9;
  
  const int steps = 100;
  const double E_step = 10.0;
  float E[steps];
  float cpp_clmb[steps];
  float fort_clmb[steps];
  float relativeDifference[steps];
  float integral = 0.0;

  float fm1 = 1.0;
  float z1 = 1.0;
  float fm2 = 82.0;
  float z2 = 50.0;
  
  for(int i = 0; i < steps; i++) {
    E[i] = E_step*i + 100.0;
    // G4double G4Incl::coulombTransm(G4double E, G4double fm1, G4double z1, G4double fm2, G4double z2)
    cpp_clmb[i] = incl->coulombTransm(double(E[i]), double(fm1), double(z1), double(fm2), double(z2));
    coulomb_transm__(&E[i], &fm1, &z1, &fm2, &z2, &fort_clmb[i]);

    if(fort_clmb[i] != 0.0) {
      relativeDifference[i] = 100.0*(cpp_clmb[i] - fort_clmb[i])/fort_clmb[i];
    }
    else {
      relativeDifference[i] = 0.0;
    }
    if(relativeDifference[i] > 2.0) {
      cout <<"TestCoulombTransm: warning! " << fort_clmb[i] << " \t \t \t " << cpp_clmb[i] << endl;
    }
    integral = integral + fabs(cpp_clmb[i] - fort_clmb[i]);
  }
  
  TCanvas *c1 = new TCanvas();

  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(steps, E, cpp_clmb);
  TGraph *fort_graph = new TGraph(steps, E, fort_clmb);

  fort_graph->GetXaxis()->SetTitle("E");
  fort_graph->GetYaxis()->SetTitle("coulomb_transm(E)");
  fort_graph->SetTitle("Function coulomb_transm");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(steps, E, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("E");
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
