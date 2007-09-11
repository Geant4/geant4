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
// $Id: TestInitIncl.cc,v 1.2 2007-09-11 13:28:42 miheikki Exp $ 
// Translation of INCL4.2/ABLA V3 
// Pekka Kaitaniemi, HIP (translation)
// Christelle Schmidt, IPNL (fission code)
// Alain Boudard, CEA (contact person INCL/ABLA)
// Aatos Heikkinen, HIP (project coordination)

// Author: Pekka Kaitaniemi <mailto:pekka.kaitaniemi@helsinki.fi>

#include "TestInitIncl.hh"

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

ClassImp(TestInitIncl)

  ///////////////////////////////////////////////////////////////////////////////////
  /* BEGIN_HTML
     <h1>Test G4Incl::initIncl</h1>

     <p>
     Test method G4Incl::initIncl
     </p>
     <p>
     <h3>Comparison between C++ and FORTRAN implementations</h3>

     <a href="logs/TestInitIncl.txt">Log file</a>
     <img src="plots/TestInitIncl.png"/>
     </p>
     END_HTML */
  
TestInitIncl::TestInitIncl()
{
  setUnitName("G4Incl::initIncl");
  setOriginalUnitName("init_incl");
  setPlotFileName("htmldoc/plots/TestInitIncl.png");
  setLogFileName("htmldoc/logs/TestInitIncl.txt");
  setLinesOfCode(8);
  setTestStatus(false);
}

TestInitIncl::~TestInitIncl()
{

}

void TestInitIncl::runMe()
{
  // This handles the actual testing.

  G4Incl *incl = new G4Incl();

  G4Calincl *calincl = (G4Calincl*) malloc(sizeof(G4Calincl));
  G4Hazard *hazard = (G4Hazard*) malloc(sizeof(G4Hazard));
  G4Spl2 *spl2 = (G4Spl2*) malloc(sizeof(G4Spl2));
  G4Dton *dton = (G4Dton*) malloc(sizeof(G4Dton));
  G4Ws *ws = (G4Ws*) malloc(sizeof(G4Ws));
  G4Mat *mat = (G4Mat*) malloc(sizeof(G4Mat));
  G4Saxw *saxw = (G4Saxw*) malloc(sizeof(G4Saxw));
  G4LightNuc *lightnuc = (G4LightNuc*) malloc(sizeof(G4LightNuc));
  G4LightGausNuc *lightgausnuc = (G4LightGausNuc*) malloc(sizeof(G4LightGausNuc));

  incl->setCalinclData(calincl);
  incl->setHazardData(hazard);
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
  
  calincl->f[0] = 207.0;
  calincl->f[1] = 82.0;
  calincl->f[3] = 1000.0;
  calincl->f[4] = 0.0;
  calincl->f[5] = 1.0;
  calincl->f[6] = 1.0;
  calincl->f[7] = 0.0;
  calincl->f[8] = 0.0;

  gCalincl->f[0] = 207.0;
  gCalincl->f[1] = 82.0;
  gCalincl->f[3] = 1000.0;
  gCalincl->f[4] = 0.0;
  gCalincl->f[5] = 1.0;
  gCalincl->f[6] = 1.0;
  gCalincl->f[7] = 0.0;
  gCalincl->f[8] = 0.0;

  hazard->ial = 39977;
  gHazard->ial = 39977;

  ws->nosurf = -1;
  gWs->nosurf = -1;

  mat->nbmat = 1;
  mat->amat[0] = 207;
  mat->zmat[0] = 82;
  gMat->nbmat = 1;
  gMat->amat[0] = 207;
  gMat->zmat[0] = 82;
  
  incl->initIncl(true);
  int ig = 1;
  init_incl__(&ig);
//   incl->initMaterial(iamat, izmat, imat);
//   init_mat__(&iamat, &izmat, &imatf);

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

  logFile << "G4Spl--------------------" << endl;
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

  setTestStatus(true);
}
