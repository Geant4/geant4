#ifndef G4GammaAnnCrossSection_h
#define G4GammaAnnCrossSection_h

#include "globals.hh"
#include "g4std/vector"
#include "G4VAnnihilationCrossSection.hh"
#include "G4ASCCrossSection.hh"

class G4GammaAnnCrossSection : public G4VAnnihilationCrossSection
{
  public:
    G4GammaAnnCrossSection();
    G4bool InCharge(G4int aCode, G4int bCode);
    G4double GetXsec(G4double s);
    
  private:
  
    G4std::vector<G4ASCCrossSection*> theGammaNucXSections;
};

#endif
