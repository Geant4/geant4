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
// $Id: TestInitMat.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestInitMat.hh"

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

ClassImp(TestInitMat)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::initMaterial</h1>

     <p>
     Test method G4Incl::initMaterial
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <a href="logs/TestInitMat.txt">Log file</a>
     <img src="plots/TestInitMat.png"/>
     </p>
     END_HTML */
  
TestInitMat::TestInitMat()
{
  setUnitName("G4Incl::initMaterial");
  setOriginalUnitName("init_mat");
  setPlotFileName("htmldoc/plots/TestInitMat.png");
  setLogFileName("htmldoc/logs/TestInitMat.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestInitMat::~TestInitMat()
{

}

void TestInitMat::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();

  G4Spl2 *spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));
  G4Dton *dton = (G4Dton*) malloc(sizeof(G4Dton));
  G4Ws *ws = (G4Ws*) malloc(sizeof(G4Ws));
  G4Mat *mat = (G4Mat*) malloc(sizeof(G4Mat));
  G4Saxw *saxw = (G4Saxw*) malloc(sizeof(G4Saxw));
  G4LightNuc *lightnuc = (G4LightNuc*) malloc(sizeof(G4LightNuc));
  G4LightGausNuc *lightgausnuc = (G4LightGausNuc*) malloc(sizeof(G4LightGausNuc));
  
  incl->setSpl2Data(spl2);
  incl->setDtonData(dton);
  incl->setWsData(ws);
  incl->setMatData(mat);
  incl->setSaxwData(saxw);
  incl->setLightNucData(lightnuc);
  incl->setLightGausNucData(lightgausnuc);

  int iamat = 207;
  int izmat = 82;
  int imat = 0;
  int imatf = 1;
  
  incl->initMaterial(iamat, izmat, imat);
  init_mat__(&iamat, &izmat, &imatf);

  ofstream logFile(getLogFileName());

  logFile << "G4Saxw--------------------" << endl;
  logFile << "Variable saxw->x[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 30; i++) {
    logFile << i << " \t \t " << gSaxw->x[0][i] << " \t \t " << saxw->x[i][0] << endl;
  }

  logFile << "Variable saxw->y[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 30; i++) {
    logFile << i << " \t \t " << gSaxw->y[0][i] << " \t \t " << saxw->y[i][0] << endl;
  }

  logFile << "Variable saxw->s[]" << endl;  
  logFile <<"i: \t \t FORTRAN: \t \t C++: " << endl;
  for(int i = 0; i < 30; i++) {
    logFile << i << " \t \t " << gSaxw->s[0][i] << " \t \t " << saxw->s[i][0] << endl;
  }

  float indices[30];
  float relativeDiff[30];
  double integral;
  double errorMarginal = 1e-3;
  
  for(int i = 0; i < 30; i++) {
    indices[i] = float(i);
    if(gSaxw->y[0][i] != 0.0) { 
      relativeDiff[i] = 100.0*((saxw->y[i][0] - gSaxw->y[0][i])/gSaxw->y[0][i]);
    }
    else {
      relativeDiff[i] = 0.0;
    }
    //    integral = integral + (spl2->a[i] + spl2->b[i] + spl2->c[i]) - (gSpl2->a[i] + gSpl2->b[i] + gSpl2->c[i]);
    integral = integral + fabs(saxw->y[i][0] - gSaxw->y[0][i]);
  }

  TGraph *diffGraph = new TGraph(30, indices, relativeDiff);
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

