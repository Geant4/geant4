#ifndef VANAPlot_h
#define VANAPlot_h

#include <iostream>
#include "Analysis/src/Particle.h"
#include "globals.hh"

class VANAPlot
{
  public:
  
    virtual G4bool Insert(G4int aPdg, G4double anX, G4double anEntryToAccumulate) = 0;
    virtual void DumpInfo(ostream &, G4String aPreFix) = 0;
    virtual void SetNevents(G4int aNumber) = 0;
    virtual G4bool Filter(ANAParticle * aParticle) = 0;
    
  private:
};

#endif
