#include "G4NeutronIsoProdCrossSections.hh"

void G4NeutronIsoProdCrossSections::
Init(ifstream & aDataSet)
{
  G4int aNumberOfPoints;
  aDataSet>>aNumberOfPoints;
  theProductionCrossSections.Init(aDataSet, aNumberOfPoints, eV);
}

G4double G4NeutronIsoProdCrossSections::
GetCrossSection(G4double anEnergy)
{
  G4double result;
  result = theProductionCrossSections.GetY(anEnergy);
  resurn result;
}

G4String G4NeutronIsoProdCrossSections::
GetProductIsotope()
{
  return theProductName;
}
