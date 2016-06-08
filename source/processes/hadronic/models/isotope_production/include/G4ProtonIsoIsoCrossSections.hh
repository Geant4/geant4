#ifndef G4ProtonIsoIsoCrossSections_h
#define G4ProtonIsoIsoCrossSections_h

#include "G4NeutronHPVector.hh"
#include "G4IsoProdCrossSections.hh"
#include "G4NeutronHPNames.hh"

class G4ProtonIsoIsoCrossSections
{

public:

  G4ProtonIsoIsoCrossSections();
  ~G4ProtonIsoIsoCrossSections();
  void Init(G4int A, G4int Z, G4double frac);
  G4double GetCrossSection(G4double anEnergy);
  G4String GetProductIsotope(G4double anEnergy);
  
  G4int GetZ() { return theZ; }
  G4int GetA() { return theA; }
  
private:
  
  G4int theNumberOfProducts;
  G4IsoProdCrossSections ** theProductionData;
  G4NeutronHPVector theCrossSection;
  G4NeutronHPNames theNames;
  G4bool hasData;
  G4int theA;
  G4int theZ;
};

#endif
