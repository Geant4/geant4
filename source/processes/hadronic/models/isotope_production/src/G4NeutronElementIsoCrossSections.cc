#include "G4NeutronElementIsoCrossSections.hh"

G4NeutronElementIsoCrossSections::
G4NeutronElementIsoCrossSections()
{
  nIsotopes = 0;
  theData = NULL;
}

G4NeutronElementIsoCrossSections::
~G4NeutronElementIsoCrossSections()
{
  for(G4int i=0; i<nIsotopes; i++)
  {
    delete theData[i];
  }
  delete theData;
}

void G4NeutronElementIsoCrossSections::
Init(const G4Element * anElement)
{
  G4int Z = anElement->GetZ();
  nIsotopes = anElement->GetNumberOfIsotopes();
  G4bool useIsotopesFromElement = true;
  if( nIsotopes == 0 ) 
  {
    nIsotopes += theStableOnes.GetNumberOfIsotopes(Z);
    useIsotopesFromElement = false;
  }
  theData = new G4NeutronIsoIsoCrossSections * [nIsotopes];
  if(useIsotopesFromElement)
  {
    for (G4int i=0; i<nIsotopes; i++)
    {
      G4int A = anElement->GetIsotope(i)->GetN();
      G4double frac = anElement->GetRelativeAbundanceVector()[i]/perCent;
      theData[i] = new G4NeutronIsoIsoCrossSections;
      theData[i]->Init(A, Z, frac);
    }
  }
  else 
  {
    G4int first = theStableOnes.GetFirstIsotope(Z);
    for(G4int i=0; i<theStableOnes.GetNumberOfIsotopes(Z); i++)
    {
      G4int A = theStableOnes.GetIsotopeNucleonCount(first+i);
      G4double frac = theStableOnes.GetAbundance(first+i);
      theData[i] = new G4NeutronIsoIsoCrossSections;
      theData[i]->Init(A, Z, frac);
    }
  }
}

G4double G4NeutronElementIsoCrossSections::
GetCrossSection(G4double anEnergy)
{
  G4double result = 0;
  for(G4int i=0; i<nIsotopes; i++)
  {
    result += theData[i]->GetCrossSection(anEnergy);
  }
  crossSectionBuffer = result;
  return result;
}

G4IsoResult * G4NeutronElementIsoCrossSections::
GetProductIsotope(G4double anEnergy)
{
  G4double running = 0;
  G4int index;
  G4double random = G4UniformRand();
  for(G4int i=0; i<nIsotopes; i++)
  {
    running += theData[i]->GetCrossSection(anEnergy);
    index = i;
    if(running/crossSectionBuffer > random) break;
  }
  G4String result = theData[index]->GetProductIsotope(anEnergy);
  G4Nucleus nucleus(theData[index]->GetA(), theData[index]->GetZ());
  G4IsoResult * theResult = new G4IsoResult(result, nucleus);
  return theResult;
}


