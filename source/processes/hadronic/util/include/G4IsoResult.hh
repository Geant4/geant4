#ifndef G4IsoResult_h
#define G4IsoResult_h 

#include "globals.hh"
#include "G4Nucleus.hh"

class G4IsoResult
{

public:

  G4IsoResult(const G4String & anIso, const G4Nucleus & aMother)
  {
    theIsoName = anIso;
    theMotherNucleus = aMother;
  }
  
  G4String GetIsotope() { return theIsoName; }
  G4Nucleus GetMotherNucleus() { return theMotherNucleus; }
  
private:
  G4IsoResult() {}
  G4IsoResult & operator = (const G4IsoResult & aResult) { return *this;}

private:
  G4String theIsoName;
  G4Nucleus theMotherNucleus;

};

#endif
