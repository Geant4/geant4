#ifndef TVANAPlot_h
#define TVANAPlot_h

#include "Analysis/src/VPlot.h"
#include "globals.hh"

template <class DataPointType>
class TVANAPlot : public VANAPlot
{
  public: 
  
  virtual G4bool Insert(G4int aPdg, G4double anX, G4double anEntryToAccumulate) = 0;
  virtual void DumpInfo(ostream &, G4String aPreFix) = 0;
};

#endif
