#ifndef G4ElementIsoCrossSections_h
#define G4ElementIsoCrossSections_h

#include "G4StableIsotopes.hh"
#include "G4IsoResult.hh"
#include "Randomize.hh"

template <class IsoIsoCrossSectionType>
class G4ElementIsoCrossSections
{

public:
  
  G4ElementIsoCrossSections();
  ~G4ElementIsoCrossSections();
  void Init(const G4Element * anElement);
  
  G4double GetCrossSection(G4double anEnergy);
  G4IsoResult * GetProductIsotope(G4double anEnergy);

private:
  
  IsoIsoCrossSectionType ** theData;
  G4int nIsotopes;
  G4StableIsotopes theStableOnes;
  G4double crossSectionBuffer;

};

#include "G4ElementIsoCrossSections.icc"

#endif
