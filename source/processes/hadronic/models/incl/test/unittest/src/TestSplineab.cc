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
// $Id: TestSplineab.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestSplineab.hh"

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

ClassImp(TestSplineab)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test splineab</h1>

     <p>
     This class tests method G4Incl::splineab.
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestSplineab.png"/>
     </p>
     END_HTML */
  
TestSplineab::TestSplineab()
{
  setUnitName("G4Incl::splineab");
  setOriginalUnitName("splineab");
  setPlotFileName("htmldoc/plots/TestSplineab.png");
  setLogFileName("htmldoc/logs/TestSplineab.log");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestSplineab::~TestSplineab()
{

}

void TestSplineab::runMe()
{
  // This handles the actual testing.
	
  G4Incl *incl = new G4Incl();

  G4Spl2 *spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));

  double errorMarginal = 1e-9;
  double integral = 0.0;
  
  const int points = 20;
  float start_point = 0.15;
  float x[points];
  Float_t cpp_y[points];
  Float_t fort_y[points];
  Float_t relativeDifference[points];
  
  float r_step = 0.1;
  int steps = 80;

  gSpl2->n = steps;
  spl2->n = steps;

  for(int i = 0; i < 100; i++) {
    spl2->x[i] = 0.0;
    spl2->y[i] = 0.0;
    spl2->a[i] = 0.0;
    spl2->b[i] = 0.0;
    spl2->c[i] = 0.0;
  }
  
  for(int i = 0; i < steps; i++) {
    gSpl2->x[i] = i*r_step;
    gSpl2->y[i] = i*r_step;
    spl2->x[i] = i*r_step;
    spl2->y[i] = i*r_step;
  }

  incl->setSpl2Data(spl2);

  // Run the functions
  incl->spl2ab();
  //  spl2ab_();

  // Make sure we have the same starting values in both FORTRAN and
  // C++ code.
  for(int i = 0; i < 100; i++) {
    gSpl2->a[i] = spl2->a[i];
    gSpl2->b[i] = spl2->b[i];
    gSpl2->c[i] = spl2->c[i];
  }

  for(int i = 0; i < points; i++) {
    x[i] = i*start_point;
    cpp_y[i] = incl->splineab(x[i]);
    fort_y[i] = splineab_(&(x[i]));
    if(fort_y[i] != 0.0) {
      relativeDifference[i] = 100.0*(cpp_y[i] - fort_y[i])/fort_y[i];
    }
    else {
      relativeDifference[i] = 0.0;
    }
    cout <<"RelDiff: " << relativeDifference[i] << endl;
    integral = integral + fabs(cpp_y[i] - fort_y[i]);
  }
  
//   ofstream logFile(getLogFileName());

//   logFile << "Variable spl2->a[]" << endl;  
//   logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
//   for(int i = 0; i < 100; i++) {
//     logFile << i << " \t \t " << gSpl2->a[i] << " \t \t " << spl2->a[i] << endl;
//   }

//   logFile << "Variable spl2->b[]" << endl;
//   cout <<"i: \t \t FORTRAN: \t \t C++: " << endl;
//   for(int i = 0; i < 100; i++) {
//     logFile << i << " \t \t " << gSpl2->b[i] << " \t \t " << spl2->b[i] << endl;
//   }

//   logFile << "Variable spl2->c[]" << endl;
//   logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
//   for(int i = 0; i < 100; i++) {
//     logFile << i << " \t \t " << gSpl2->c[i] << " \t \t " << spl2->c[i] << endl;
//   }
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);

  TGraph *cpp_graph = new TGraph(points, x, cpp_y);
  TGraph *fort_graph = new TGraph(points, x, fort_y);

  fort_graph->GetXaxis()->SetTitle("x");
  fort_graph->GetYaxis()->SetTitle("splineab(x)");
  fort_graph->SetTitle("Function splineab");
  fort_graph->Draw("al");
  cpp_graph->Draw("p, same");	

  c1->cd(4);
  TGraph *diff_graph = new TGraph(points, x, relativeDifference);
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

  if(integral < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
