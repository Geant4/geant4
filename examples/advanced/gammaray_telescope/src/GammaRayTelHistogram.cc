// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelHistogram.cc,v 1.1 2000-11-15 20:27:41 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      AIDA interface
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelHistogram  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#ifdef G4HIS_USE_AID
#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelHistogram.hh"

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

GammaRayTelHistogram::GammaRayTelHistogram(GammaRayTelDetectorConstruction* GammaRayTelDC): 
GammaRayTelDetector(GammaRayTelDC), hFact(0), vFact(0), pl1(0), pl2(0) {
  
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

  pl1 = createIPlotter(); // create plotter with 1 zone
  if (pl1 == 0) {
    cerr << "ERROR: could not create plotter !" << endl;
    exit(-1);
  }

  pl2 = createIPlotter(); // create plotter with 1 zone
  if (pl2 == 0) {
    cerr << "ERROR: could not create plotter !" << endl;
    exit(-1);
  }

  book();
}

GammaRayTelHistogram::~GammaRayTelHistogram() {
  delete hFact;
  delete vFact;
  delete pl1;
  //  delete pl2;
}

bool GammaRayTelHistogram::book() {
  float sizexy, sizez;
  
  sizexy = GammaRayTelDetector->GetTKRSizeXY();
  sizez = GammaRayTelDetector->GetTKRSizeZ();

  // Create 1D histogram for the kinetic Energy
  // No need to specify min/max for dynamical histograms. 
  // The data of the "fill"s are stored (up to 10^6) as they come in.
  // When reaching 10^6 "fill"s the histogram is "frozen" and behaves "normal"
  IHistogram1D *dynHist = hFact->dynamic1D(10,"kinetic Energy",100);
  if (dynHist == 0) {
    cerr << "ERROR: could not create dynamical 1D histogram for kinetic energy!" << endl;
    return false;
  }
  h1dList.push_back(dynHist);
  
  // Book the histogram for the 2D position. 
  // Instead of using a scatter plot, just book enough bins ...
  IHistogram2D *pos1 = hFact->create2D(100,"Tracker Hits XZ", 
				       100, -sizez/2, sizez/2, 
				       100, -sizexy/2, sizexy/2);
  h2dList.push_back(pos1);
  IHistogram2D *pos2 = hFact->create2D(200,"Tracker Hits YZ", 
				       100, -sizez/2, sizez/2, 
				       100, -sizexy/2, sizexy/2);
  h2dList.push_back(pos2);
  return true;
}

bool GammaRayTelHistogram::finish() {
  
  pl1->zone(2,1);		// create 2 zones in page

  plot(h1dList[0]);
  
  pl1->refresh();

  pl1->psPrint("g4-mini.ps");

  return true;
}

bool GammaRayTelHistogram::plot(IHistogram1D * hist) {

  // create a VectorOfPoints from histogram :
  IVector * v = vFact->from1D(hist);
  if (v == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl1->plot(v);
  
  //delete v;			// don't forget this !

  return true;

}

bool GammaRayTelHistogram::plot1(IHistogram2D * hist) {

  // create a VectorOfPoints from histogram :
  IVector * v = vFact->from2D(hist);
  if (v == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl1->plot(v);

  //delete v;			// don't forget this !

  return true;

}

bool GammaRayTelHistogram::plot2(IHistogram2D * hist) {

  // create a VectorOfPoints from histogram :
  IVector * v = vFact->from2D(hist);
  if (v == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl2->plot(v);

  //delete v;			// don't forget this !

  return true;
  
}
#endif











