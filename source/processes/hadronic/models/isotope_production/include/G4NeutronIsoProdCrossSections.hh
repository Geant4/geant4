#ifndef G4NeutronIsoProdCrossSections_h
#define G4NeutronIsoProdCrossSections_h

#include "globals.hh"
#include "G4NeutronHPVector.hh"

class G4NeutronIsoProdCrossSections
{
public:
  G4NeutronIsoProdCrossSections(G4String aString)
  { theProductName=aString; }
  void Init(ifstream & aDataSet);
  G4double GetProductionCrossSection(G4double anEnergy);
  G4String GetProductIsotope();

private:

  G4String theProductName;
  G4NeutronHPVector theProductionCrossSections;
};

#endif
