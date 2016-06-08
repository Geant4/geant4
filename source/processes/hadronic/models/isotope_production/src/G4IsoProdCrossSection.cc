#include "G4IsoProdCrossSections.hh"

void G4IsoProdCrossSections::
Init(G4std::ifstream & aDataSet)
{
  G4int aNumberOfPoints;
  aDataSet>>aNumberOfPoints;
  theProductionCrossSections.Init(aDataSet, aNumberOfPoints, eV);
}

G4double G4IsoProdCrossSections::
GetProductionCrossSection(G4double anEnergy)
{
  G4double result;
  result = theProductionCrossSections.GetY(anEnergy);
  return result;
}

G4String G4IsoProdCrossSections::
GetProductIsotope()
{
  return theProductName;
}
