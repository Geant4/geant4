

#include "XrayTelHistogram.hh"

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

XrayTelHistogram::XrayTelHistogram() : hFact(0), vFact(0), pl(0) {

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
  delete hFact;
  delete vFact;
  delete pl;
}

bool XrayTelHistogram::book() {

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
  IHistogram2D *pos = hFact->create2D(100,"position of hist", 200, -50, 50, 200, -50, 50);
  h2dList.push_back(pos);

  return true;
}

bool XrayTelHistogram::finish() {

  pl->zone(2,1);		// create 2 zones in page

  plot(h1dList[0]);

  pl->refresh();		// not really needed here, but that's needed if the plotting window was covered ...

  pl->psPrint("g4-mini.ps");

  return true;
}

bool XrayTelHistogram::plot(IHistogram1D * hist) {

  // create a VectorOfPoints from histogram :
  IVector * v = vFact->from1D(hist);
  if (v == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl->plot(v);

  //delete v;			// don't forget this !

  return true;

}

bool XrayTelHistogram::plot(IHistogram2D * hist) {

  // create a VectorOfPoints from histogram :
  IVector * v = vFact->from2D(hist);
  if (v == 0) {
    cerr << "ERROR: could not create vector !" << endl;
    return false;
  }
  // now plot the vector
  pl->plot(v);

  //delete v;			// don't forget this !

  return true;

}


