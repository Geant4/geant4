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
// * MODULE:            XrayTelHistogram.hh                             *     
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

#include "g4std/vector"

#include "interfaces/IPlotter.h"
#include "interfaces/IHistogram.h"
#include "interfaces/IHistogram1D.h"
#include "interfaces/IHistogram2D.h"

class G4SteppingManager;

class IHistoFactory;
class IVectorFactory;
class IPlotter;

class XrayTelHistogram {

public:
  XrayTelHistogram();
  virtual ~XrayTelHistogram();

  bool book();
  bool finish();

  bool plot(IHistogram1D *);
  bool plot(IHistogram2D *);

  IPlotter * getPlotter() { return pl; }

public:
  // callback functions from various places, distinguished by argument passed
  // see visitor pattern
  void analyze(const G4SteppingManager *);
  // ... what else ?

public:
  // to have access to the histograms 
  vector<IHistogram1D *> * getH1DList() { return &h1dList; }
  vector<IHistogram2D *> * getH2DList() { return &h2dList; }

private:
  IHistoFactory  *hFact;
  IVectorFactory *vFact;
  IPlotter       *pl;

  vector<IHistogram1D *> h1dList;
  vector<IHistogram2D *> h2dList;

  IVector * vCache;

};
