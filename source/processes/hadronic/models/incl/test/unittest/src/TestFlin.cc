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
// $Id: TestFlin.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestFlin.hh"

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

ClassImp(TestFlin)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test flin</h1>

     <p>
     This class tests method G4Incl::interpolateFunction
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestFlin.png"/>
     </p>
     <p>
     Extra logs:
     <a href="logs/TestFlinExtraLogCpp.txt">C++</a> <a href="logs/TestFlinExtraLogFortran.txt">FORTRAN</a>
     </p>
     END_HTML */
  
TestFlin::TestFlin()
{
  setUnitName("G4Incl::interpolateFunction");
  setOriginalUnitName("flin");
  setPlotFileName("htmldoc/plots/TestFlin.png");
  setLogFileName("htmldoc/logs/TestFlin.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestFlin::~TestFlin()
{

}

void TestFlin::runMe()
{
  // This handles the actual testing.
	
  // Create a new G4Incl instance.
  G4Incl *incl = new G4Incl();
  
  // Set starting values:
  Float_t r0 = 0.1;
  Float_t adif = 1.0;
  Float_t integral = 0.0;

  // Set error marginal
  const Float_t errorMarginal = 1e-9;

  // Set up arrays for data collection:
  const Int_t points = 10;
  const Float_t r_step = 1.0;
  Float_t r[points];
  r[0] = 0.0;
  Float_t cpp_f_r[points];
  Float_t fort_f_r[points];
  Float_t relativeDifference[points];

  // Starting values:
  G4Saxw *saxw = (G4Saxw*) malloc(sizeof(G4Saxw));
  incl->setSaxwData(saxw);
  int flinpoint;

  saxw->k = 0;
  gSaxw->k = 1;
  saxw->n = 30;
  gSaxw->n = 30;
  for(Int_t i = 0; i < 30; i++) {
    saxw->x[i][0] = i*r_step + 0.05;
    saxw->y[i][0] = 2*(i*r_step) + (i*r_step)*(i*r_step);
    //saxw->y[i][0] = 2*(i*r_step);
    gSaxw->x[0][i] = saxw->x[i][0];
    gSaxw->y[0][i] = saxw->y[i][0];
  }

  int material = 1;
  //  flin2_(&material);
  incl->firstDerivative(0);
  
  for(Int_t i = 0; i < gSaxw->n; i++) {
    gSaxw->s[0][i] = saxw->s[i][0];
  }
  
  ofstream cppnumbers("htmldoc/logs/TestFlinExtraLogCpp.txt");
  ofstream fortnumbers("htmldoc/logs/TestFlinExtraLogFortran.txt");
  for(Int_t k = 0; k < 500; k++) {
    cppnumbers << "Material (k) = " << k << endl;
    fortnumbers << "Material (k) = " << k << endl;
    cppnumbers << "x[n][k] \t \t \t y[n][k]" << endl;
    fortnumbers << "x[n][k] \t \t \t y[n][k]" << endl;
    cppnumbers.precision(3);
    fortnumbers.precision(3);
    for(Int_t n = 0; n < 30; n++) {
      cppnumbers << saxw->x[n][k]  << "\t \t \t" << saxw->y[n][k] << " \t \t \t" << saxw->s[n][k] << endl;
      fortnumbers << gSaxw->x[k][n]  << "\t \t \t" << gSaxw->y[k][n] << " \t \t \t" << gSaxw->s[k][n] << endl;
    }
    cppnumbers << endl;
    fortnumbers << endl;
  }
  cppnumbers.close();
  fortnumbers.close();

  for(Int_t i = 1; i < points; i++) {
    r[i] = r[i-1] + r_step;	
    cpp_f_r[i] = incl->interpolateFunction(r[i]);
    fort_f_r[i] = flin_(&r[i]);
    
    relativeDifference[i] = 100.0*(cpp_f_r[i] - fort_f_r[i])/fort_f_r[i];
    integral = integral + abs(cpp_f_r[i] - fort_f_r[i]);
  }

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, r, cpp_f_r);
  TGraph *fort_graph = new TGraph(points, r, fort_f_r);

  fort_graph->GetXaxis()->SetTitle("x");
  fort_graph->GetYaxis()->SetTitle("flin(x)");
  fort_graph->SetTitle("Function flin");

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
