#ifndef G4NeutronIsoIsoCrossSections_h
#define G4NeutronIsoIsoCrossSections_h

#include "G4NeutronHPVector.hh"
#include "G4NeutronIsoProdCrossSections.hh"
#include "G4NeutronHPNames.hh"

class G4NeutronIsoIsoCrossSections
{

public:

  G4NeutronIsoIsoCrossSections();
  ~G4NeutronIsoIsoCrossSections();
  void Init(G4int A, G4int Z, G4double frac);
  G4double GetCrossSection(G4double anEnergy);
  G4String GetProductIsotope(G4double anEnergy);
  
  G4int GetZ() { return theZ; }
  G4int GetA() { return theA; }
  
private:
  
  G4int theNumberOfProducts;
  G4NeutronIsoProdCrossSections ** theProductionData;
  G4NeutronHPVector theCrossSection;
  G4NeutronHPVector theInelasticCrossSection;
  G4NeutronHPVector theNNprimeCrossSection;
  G4NeutronHPNames theNames;
  G4bool hasData;
  G4int theA;
  G4int theZ;
};

#endif
