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
// $Id: TestClmb1.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestClmb1.hh"

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

ClassImp(TestClmb1)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::clmb1</h1>

     <p>
     Test for G4Incl::clmb1
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestClmb1.png"/>
     </p>
     END_HTML */
  
TestClmb1::TestClmb1()
{
  setUnitName("G4Incl::clmb1");
  setOriginalUnitName("clmb1");
  setPlotFileName("htmldoc/plots/TestClmb1.png");
  setLogFileName("htmldoc/logs/TestClmb1.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestClmb1::~TestClmb1()
{

}

void TestClmb1::runMe()
{
  // This handles the actual testing.

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
  
  double eta = 40.0;
  double t1c = 1.0;
  double t1fort = 1.0;
  
  for(int i = 0; i < steps; i++) {
    rho[i] = rho_step*i + 1.0;
    cout <<"Calling C++" << endl;
    cpp_clmb[i] = incl->clmb1(rho[i], eta, &t1c);
    cout <<"Calling fortran" << endl;
    cout <<"rho[i] : " << rho[i] << " eta: " << eta << " t1: " << t1fort << endl;
    fort_clmb[i] = clmb1_(&rho[i], &eta, &t1fort);
    cout <<"Fortran call complete. Result: " << fort_clmb[i] << endl;

    if(fort_clmb[i] != 0.0) {
      relativeDifference[i] = 100.0*(cpp_clmb[i] - fort_clmb[i])/fort_clmb[i];
    }
    else {
      relativeDifference[i] = 0.0;
    }
    if(relativeDifference[i] > 2.0) {
      cout <<"TestClmb1: Relative diff. warning: " << relativeDifference[i] << " rho: " << rho[i] <<" eta: " << eta << endl;
    }
    integral = integral + fabs(cpp_clmb[i] - fort_clmb[i]);

    cout <<fort_clmb[i] << " \t \t \t " << cpp_clmb[i] << endl;
  }
  
  TCanvas *c1 = new TCanvas();

  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(steps, rho, cpp_clmb);
  TGraph *fort_graph = new TGraph(steps, rho, fort_clmb);

  cout <<"Graphs generated..." << endl;
  
  fort_graph->GetXaxis()->SetTitle("#rho");
  fort_graph->GetYaxis()->SetTitle("clmb1(#rho)");
  fort_graph->SetTitle("Function clmb1");

  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  cout <<"Graph 1 ready." << endl;

  c1->cd(4);
  TGraph *diff_graph = new TGraph(steps, rho, relativeDifference);
  diff_graph->SetTitle("Relative difference C++/FORTRAN");
  diff_graph->GetXaxis()->SetTitle("#rho");
  diff_graph->GetYaxis()->SetTitle("Relative difference (%)");

  diff_graph->Draw("ap");

  cout <<"Graph 2 ready." << endl;
  
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
