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
// $Id: TestXabs2.cc,v 1.2 2007-09-11 13:28:43 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestXabs2.hh"

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

ClassImp(TestXabs2)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::xabs2</h1>

     <p>
     Test method G4Incl::xabs2.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestXabs2.png"/>
     </p>
     END_HTML */
  
TestXabs2::TestXabs2()
{
  setUnitName("G4Incl::xabs2");
  setOriginalUnitName("xabs2");
  setPlotFileName("htmldoc/plots/TestXabs2.png");
  setLogFileName("htmldoc/logs/TestXabs2.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestXabs2::~TestXabs2()
{

}

void TestXabs2::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();


  const double errorMarginal = 1e-2;
  
  const int steps = 100;
  const double E_step = 10.0;
  double E[steps];
  double cpp_result[steps];
  double fort_result[steps];
  float temp;
  double relativeDifference[steps];
  double integral = 0.0;

  double zp = 1.0;
  double ap = 1.0;
  double zt = 50.0;
  double at = 82.0;

    //  G4double G4Incl::xabs2(G4double zp_p, G4double ap_p, G4double zt_p, G4double at_p, G4double ep_p)
  //void xabs2_(double *zp, double *ap, double *zt, double *at, double *ep, float *sig);

  for(int i = 0; i < steps; i++) {
    E[i] = E_step*i + 100.0;
    // G4double G4Incl::coulombTransm(G4double E, G4double fm1, G4double z1, G4double fm2, G4double z2)
    cpp_result[i] = incl->xabs2(zp, ap, zt, at, E[i]);
    xabs2_(&zp, &ap, &zt, &at, &E[i], &temp);
    fort_result[i] = temp;
    
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
  fort_graph->GetYaxis()->SetTitle("xabs2(E)");
  fort_graph->SetTitle("Function xabs2");

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
