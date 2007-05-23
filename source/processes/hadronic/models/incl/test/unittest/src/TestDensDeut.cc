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
// $Id: TestDensDeut.cc,v 1.1 2007-05-23 09:56:26 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestDensDeut.hh"

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

ClassImp(TestDensDeut)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::densDeut</h1>

     <p>
     Test method G4Incl::densDeut
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <img src="plots/TestDensDeut.png"/>
     </p>
     END_HTML */
  
TestDensDeut::TestDensDeut()
{
  setUnitName("G4Incl::densDeut");
  setOriginalUnitName("dens_deut");
  setPlotFileName("htmldoc/plots/TestDensDeut.png");
  setLogFileName("htmldoc/logs/TestDensDeut.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestDensDeut::~TestDensDeut()
{

}

void TestDensDeut::runMe()
{
  // This handles the actual testing.
	
  G4Incl *incl = new G4Incl();

  G4Spl2 *spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));
  G4Dton *dton = (G4Dton*) malloc(sizeof(G4Dton));
  G4Ws *ws = (G4Ws*) malloc(sizeof(G4Ws));
  
  incl->setSpl2Data(spl2);
  incl->setDtonData(dton);
  incl->setWsData(ws);
  
  incl->densDeut();
  dens_deut__();

  const double errorMarginal = 1e-3;
  float indices[100];
  float relativeDiff[100];
  double integral;
  
  ofstream logFile(getLogFileName());

  logFile << "Variable spl2->x[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 100; i++) {
    logFile << i << " \t \t " << gSpl2->x[i] << " \t \t " << spl2->x[i] << endl;
  }

  logFile << "Variable spl2->y[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 100; i++) {
    logFile << i << " \t \t " << gSpl2->y[i] << " \t \t " << spl2->y[i] << endl;
  }

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
    if((gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]) != 0.0) { 
      relativeDiff[i] = 100.0*((spl2->a[i] + spl2->b[i] + spl2->c[i]) - (gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]))/(gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]);
    }
    else {
      relativeDiff[i] = 0.0;
    }
    //    integral = integral + (spl2->a[i] + spl2->b[i] + spl2->c[i]) - (gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]);
    integral = integral + fabs(spl2->x[i] - gSpl2->x[i]);
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

  cout <<"Integral: " << integral << endl;
  
  if(integral < errorMarginal) {
    setTestStatus(true);
  }
  else {
    setTestStatus(false);
  }
}
