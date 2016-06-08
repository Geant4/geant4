// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VAnalysisSystem.hh,v 1.8 2000/11/16 13:44:46 barrand Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
// 

#ifndef G4VANALYSISSYSTEM_HH
#define G4VANALYSISSYSTEM_HH

#ifdef G4ANALYSIS_BUILD

#include "g4rw/cstring.h"

class IHistogram;
class IHistogramFactory;

#ifdef G4ANALYSIS_NON_AIDA
class ICloudFactory;
class ICloud;
class ITuple;
#endif

class G4VAnalysisSystem {
public:
  virtual ~G4VAnalysisSystem() {}
  virtual const G4String& GetName() const = 0;
  virtual IHistogramFactory* GetHistogramFactory() = 0;
  //
  virtual void Store(IHistogram* = 0,const G4String& = "") = 0;
  virtual void Plot(IHistogram*) = 0;
  //
#ifdef G4ANALYSIS_NON_AIDA
  virtual ICloudFactory* GetCloudFactory() = 0;
  virtual void Plot(ICloud*) = 0;
  virtual ITuple* CreateTuple(const G4String&,const G4String&) = 0;
#endif
};

#endif

#endif
