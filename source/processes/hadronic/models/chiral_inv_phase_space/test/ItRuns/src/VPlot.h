#ifndef VANAPlot_h
#define VANAPlot_h

#include <iostream>

class VANAPlot
{
  public:
  
  virtual G4bool Insert(G4int aPdg, G4double anX, G4double anEntryToAccumulate) = 0;
  virtual void DumpInfo(ostream &) = 0;
};

#endif
