// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelHistogram.hh,v 1.1 2000-11-15 20:27:39 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//     
//      AIDA interface
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelHistogram  ------
//           by R.Giannitrapani, F. Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#ifdef G4HIS_USE_AIDA
#ifndef GammaRayTelHistogram_h
#define GammaRayTelHistogram_h 1
#include <vector>

#include "interfaces/IPlotter.h"

#include "interfaces/IHistogram.h"
#include "interfaces/IHistogram1D.h"
#include "interfaces/IHistogram2D.h"

class GammaRayTelDetectorConstruction;
class IHistoFactory;
class IVectorFactory;
class IPlotter;

class GammaRayTelHistogram 
{
public:
  GammaRayTelHistogram(GammaRayTelDetectorConstruction*);
  virtual ~GammaRayTelHistogram();
  
  bool book();
  bool finish();

  bool plot(IHistogram1D *);
  bool plot1(IHistogram2D *);
  bool plot2(IHistogram2D *);

  IPlotter * getPlotter1() { return pl1; }
  IPlotter * getPlotter2() { return pl2; }

public:
  // to have access to the histograms 
  vector<IHistogram1D *> * getH1DList() { return &h1dList; }
  vector<IHistogram2D *> * getH2DList() { return &h2dList; }

private:
  IHistoFactory  *hFact;
  IVectorFactory *vFact;
  IPlotter       *pl1, *pl2;

  vector<IHistogram1D *> h1dList;
  vector<IHistogram2D *> h2dList;

  GammaRayTelDetectorConstruction*    GammaRayTelDetector;
};


#endif
#endif






