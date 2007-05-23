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
// $Id: TestSpl2Ab.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestSpl2Ab.hh"

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

ClassImp(TestSpl2Ab)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::spl2ab</h1>

     <p>
     This class tests method G4Incl::spl2ab
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestSpl2Ab.png"/>
     </p>
     <p>
     A more detailed <a href="logs/TestSpl2Ab.txt">log file</a>.
     </p>
     END_HTML */
  
TestSpl2Ab::TestSpl2Ab()
{
  setUnitName("G4Incl::spl2Ab");
  setOriginalUnitName("spl2ab");
  setPlotFileName("htmldoc/plots/TestSpl2Ab.png");
  setLogFileName("htmldoc/logs/TestSpl2Ab.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestSpl2Ab::~TestSpl2Ab()
{

}

void TestSpl2Ab::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();

  G4Spl2 *spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));

  const Double_t errorMarginal = 1e-9;
  Double_t integral = 0.0;
  
  float r_step = 0.1;
  int steps = 80;
  Float_t indices[100];
  Float_t relativeDiff[100];

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
  spl2ab_();

  ofstream logFile(getLogFileName());

  logFile << "Variable spl2->a[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 100; i++) {
    logFile << i << " \t \t " << gSpl2->a[i] << " \t \t " << spl2->a[i] << endl;
  }

  logFile << "Variable spl2->b[]" << endl;
  cout <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 100; i++) {
    logFile << i << " \t \t " << gSpl2->b[i] << " \t \t " << spl2->b[i] << endl;
  }

  logFile << "Variable spl2->c[]" << endl;
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 100; i++) {
    logFile << i << " \t \t " << gSpl2->c[i] << " \t \t " << spl2->c[i] << endl;
  }

  for(int i = 0; i < 100; i++) {
    indices[i] = float(i);
    relativeDiff[i] = 100.0*((spl2->a[i] + spl2->b[i] + spl2->c[i]) - (gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]))/(gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]);
    integral = integral + (spl2->a[i] + spl2->b[i] + spl2->c[i]) - (gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]);
  }

  TGraph *diffGraph = new TGraph(100, indices, relativeDiff);
  diffGraph->SetTitle("Relative difference C++/FORTRAN");
  diffGraph->GetXaxis()->SetTitle("Array index");
  diffGraph->GetYaxis()->SetTitle("Relative difference (%)");

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);

  c1->cd(4);
  diffGraph->Draw("ap");
  
  c1->SaveAs(getPlotFileName());

  if(integral < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
