#ifndef G4IsoProdCrossSections_h
#define G4IsoProdCrossSections_h

#include "globals.hh"
#include "G4NeutronHPVector.hh"

class G4IsoProdCrossSections
{
public:
  G4IsoProdCrossSections(G4String aString)
  { theProductName=aString; }
  void Init(G4std::ifstream & aDataSet);
  G4double GetProductionCrossSection(G4double anEnergy);
  G4String GetProductIsotope();

private:

  G4String theProductName;
  G4NeutronHPVector theProductionCrossSections;
};

#endif
