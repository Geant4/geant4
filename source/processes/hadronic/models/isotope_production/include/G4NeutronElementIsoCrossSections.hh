#ifndef G4NeutronElementIsoCrossSections_h
#define G4NeutronElementIsoCrossSections_h

#include "G4NeutronIsoIsoCrossSections.hh"
#include "G4StableIsotopes.hh"
#include "G4IsoResult.hh"
#include "Randomize.hh"

class G4NeutronElementIsoCrossSections
{

public:
  
  G4NeutronElementIsoCrossSections();
  ~G4NeutronElementIsoCrossSections();
  void Init(const G4Element * anElement);
  
  G4double GetCrossSection(G4double anEnergy);
  G4IsoResult * GetProductIsotope(G4double anEnergy);

private:
  
  G4NeutronIsoIsoCrossSections ** theData;
  G4int nIsotopes;
  G4StableIsotopes theStableOnes;
  G4double crossSectionBuffer;

};

#endif
