// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelHistogram.cc                             *     
// * -------                                                            *
// *                                                                    *
// * Version:                                                           *
// * Date:                                                              *
// * Author:            A Pfeiffer                                      *
// * Organisation:      CERN/IT, Geneva                                 *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
//
// **********************************************************************

#include "G4ios.hh"

#include "XrayTelSteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"

//#include "XrayTelHistogram.hh"

#include <iostream>
#include <string>
#include <stdlib.h>

#include "AIDA_HTL/AIDAHistograms.h"

#include "AIDA_Plotter/AIDAPlotter.h"
// get access to the VectorOfPoints 
#include "interfaces/IVector.h"
#include "interfaces/IVectorFactory.h"
// .. and the Annotation
#include "interfaces/IAnnotation.h"
#include "interfaces/IAnnotationFactory.h"

XrayTelHistogram::XrayTelHistogram() : hFact(0), vFact(0), pl(0), vCache(0) {

  // once we inherit from G4AnalysisManager, the following will change
  // as the initialisation/loading of libs will be done there ...

  hFact = createIHistoFactory();
  if (hFact == 0) {
    cerr << "ERROR: could not create histoFactory !" << endl;
    exit(-1);
  }

  vFact = createIVectorFactory();
  if (vFact == 0) {
    cerr << "ERROR: could not create vector factory !" << endl;
    exit(-1);
  }

  pl = createIPlotter(); // create plotter with 1 zone
  if (pl == 0) {
    cerr << "ERROR: could not create plotter !" << endl;
    exit(-1);
  }

  book();
}

XrayTelHistogram::~XrayTelHistogram() {

  // ... last chance to do things ... :-)

  if (h1dList[0] != 0) {
    plot(h1dList[0]);
    pl->refresh();
    pl->psPrint("kineticEnergy.ps");
    vCache->toAscii("kineticEnergy.vec");
    G4cerr << "kinetic energy plot created" << G4endl;
  } else {
    G4cerr << "kinetic energy histo not there !?!?" << G4endl;
  }

  // ... cleaning up ...

  delete hFact;
  delete vFact;
  delete pl;

  delete vCache;
}

void XrayTelHistogram::analyze(const G4SteppingManager *pSM) {

  G4Track* fTrack = pSM->GetTrack();

  Hep3Vector pos = fTrack->GetPosition();

  IHistogram1D * kinE;
  kinE = h1dList[0]; // 1D-histo # 0 is kinetic engergy, could also be a direct pointer to a histo ...
  if (kinE == 0) {
    G4cerr << "XrayTelSteppingAction::UserSteppingAction> could not find histo for kinetic energy !" << G4endl;
  } else {
    kinE->fill(fTrack->GetKineticEnergy());
  }
  
  IHistogram2D * posHist;
  posHist = h2dList[0];// 2D-histo # 0 is pos, could also be a direct pointer to a histo ...
  if (posHist == 0) {
    G4cerr << "XrayTelSteppingAction::UserSteppingAction> could not find histo for position !" << G4endl;
  } else {
    posHist->fill(pos.y(), pos.z());
    plot(posHist);
    pl->refresh();
    pl->psPrint("posplot.ps");
  }
  
}
bool XrayTelHistogram::book() {

  // Create 1D histogram for the kinetic Energy
  // No need to specify min/max for dynamical histograms. 
  // The data of the "fill"s are stored (up to 10^6) as they come in.
  // When reaching 10^6 "fill"s the histogram is "frozen" and behaves "normal"
  IHistogram1D *dynHist = hFact->create1D(10,"kinetic Energy",100,0.,0.5);
  if (dynHist == 0) {
    cerr << "ERROR: could not create dynamical 1D histogram for kinetic energy!" << endl;
    return false;
  }
  h1dList.push_back(dynHist);

  for (int i=0; i<100; i++) {
    dynHist->fill(.3/100.*i);
  }

  // Book the histogram for the 2D position. 
  // Instead of using a scatter plot, just book enough bins ...
  IHistogram2D *pos = hFact->create2D(100,"position of hist", 200, -50, 50, 200, -50, 50);
  h2dList.push_back(pos);

  return true;
}

bool XrayTelHistogram::finish() {

  return true;
}

bool XrayTelHistogram::plot(IHistogram1D * hist) {

  delete vCache;
  // create a VectorOfPoints from histogram :
  vCache = vFact->from1D(hist);
  if (vCache == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl->plot(vCache);
  return true;

}

bool XrayTelHistogram::plot(IHistogram2D * hist) {

  delete vCache;
  // create a VectorOfPoints from histogram :
  vCache = vFact->from2D(hist);
  if (vCache== 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl->plot(vCache);

  return true;

}


