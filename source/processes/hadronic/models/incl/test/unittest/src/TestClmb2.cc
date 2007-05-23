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
// $Id: TestClmb2.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestClmb2.hh"

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

ClassImp(TestClmb2)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test Clmb2</h1>

     <p>
     Test for the method G4Incl::clmb2.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestClmb2.png"/>
     </p>
     END_HTML */
  
TestClmb2::TestClmb2()
{
  setUnitName("G4Incl::clmb2");
  setOriginalUnitName("clmb2");
  setPlotFileName("htmldoc/plots/TestClmb2.png");
  setLogFileName("htmldoc/logs/TestClmb2.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestClmb2::~TestClmb2()
{

}

void TestClmb2::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();

  const double errorMarginal = 1e-9;
  
  const int steps = 100;
  const double rho_step = 1.0;
  double rho[steps];
  double cpp_clmb[steps];
  double fort_clmb[steps];
  double relativeDifference[steps];
  double integral = 0.0;
  
  double eta = 1.0;
  double t1c = 1.0;
  double t1fort = 1.0;
  
  for(int i = 0; i < steps; i++) {
    rho[i] = rho_step*i;
    cpp_clmb[i] = incl->clmb2(rho[i], eta, &t1c);
    fort_clmb[i] = clmb2_(&rho[i], &eta, &t1fort);

    if(fort_clmb[i] != 0.0) {
      relativeDifference[i] = 100.0*(cpp_clmb[i] - fort_clmb[i])/fort_clmb[i];
    }
    else {
      relativeDifference[i] = 0.0;
    }
    if(relativeDifference[i] > 2.0) {
      cout <<"TestClmb2: Relative diff. warning: " << relativeDifference[i] << " rho: " << rho[i] <<" eta: " << eta << endl;
    }
    integral = integral + fabs(cpp_clmb[i] - fort_clmb[i]);
  }
  
  TCanvas *c1 = new TCanvas();

  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(steps, rho, cpp_clmb);
  TGraph *fort_graph = new TGraph(steps, rho, fort_clmb);

  fort_graph->GetXaxis()->SetTitle("#rho");
  fort_graph->GetYaxis()->SetTitle("clmb2(#rho)");
  fort_graph->SetTitle("Function clmb2");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(steps, rho, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("#rho");
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
