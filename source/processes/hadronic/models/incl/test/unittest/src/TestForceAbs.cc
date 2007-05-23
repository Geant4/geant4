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
// $Id: TestForceAbs.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestForceAbs.hh"

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

ClassImp(TestForceAbs)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::forceAbs</h1>

     <p>
     Test method G4Incl::forceAbs
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestForceAbs.png"/>
     </p>
     END_HTML */
  
TestForceAbs::TestForceAbs()
{
  setUnitName("G4Incl::forceAbs");
  setOriginalUnitName("force_abs");
  setPlotFileName("htmldoc/plots/TestForceAbs.png");
  setLogFileName("htmldoc/logs/TestForceAbs.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestForceAbs::~TestForceAbs()
{

}

void TestForceAbs::runMe()
{
  // This handles the actual testing.
	

  G4Incl *incl = new G4Incl();


  const double errorMarginal = 1e-2;
  
  const int steps = 100;
  const double E_step = 10.0;
  float E[steps];
  float cpp_result[steps];
  float fort_result[steps];
  float temp;
  float relativeDifference[steps];
  double integral = 0.0;

    //  G4double G4Incl::xabs2(G4double zp_p, G4double ap_p, G4double zt_p, G4double at_p, G4double ep_p)
  //void xabs2_(double *zp, double *ap, double *zt, double *at, double *ep, float *sig);

  int iprojo = 1;
  float at = 207.0;
  float zt = 82.0;
  float bmax = 1.0;
  float pt = 10.0;
  
  for(int i = 0; i < steps; i++) {
    E[i] = E_step*i + 100.0;
    cpp_result[i] = incl->forceAbs(iprojo, at, zt, E[i], bmax, pt);
    force_abs__(&iprojo, &at, &zt, &E[i], &bmax, &pt, &temp);
    fort_result[i] = temp;
    cout <<"Fortran: " << temp << " \t \t C++: " << cpp_result[i] << endl; 
    if(fort_result[i] != 0.0) {
      relativeDifference[i] = 100.0*(cpp_result[i] - fort_result[i])/fort_result[i];
    }
    else {
      relativeDifference[i] = 0.0;
    }
    if(relativeDifference[i] > 2.0) {
      cout <<"TestCoulombTransm: warning! " << fort_result[i] << " \t \t \t " << cpp_result[i] << endl;
    }
    integral = integral + fabs(cpp_result[i] - fort_result[i]);
  }
  
  TCanvas *c1 = new TCanvas();

  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(steps, E, cpp_result);
  TGraph *fort_graph = new TGraph(steps, E, fort_result);

  fort_graph->GetXaxis()->SetTitle("E");
  fort_graph->GetYaxis()->SetTitle("forceAbs(E)");
  fort_graph->SetTitle("Function forceAbs");

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

  cout <<"Integral : " << integral << endl;
  
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
