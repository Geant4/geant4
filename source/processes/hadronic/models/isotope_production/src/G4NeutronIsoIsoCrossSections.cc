#include "G4NeutronIsoIsoCrossSections.hh"
#include "G4NeutronHPDataUsed.hh"

G4NeutronIsoIsoCrossSections::
G4NeutronIsoIsoCrossSections()
: theCrossSection(), theNames()
{
  theProductionData = NULL;
  hasData = false;
  theNumberOfProducts = 0;
  theZ = 0;
  theA = 0;
}

G4NeutronIsoIsoCrossSections::
~G4NeutronIsoIsoCrossSections()
{
  for(G4int i=0; i<theNumberOfProducts; i++)
  {
    delete theProductionData[i];
  }
  delete theProductionData;
}

void G4NeutronIsoIsoCrossSections::
Init(G4int A, G4int Z, G4double frac)
{
  // First transmution scattering cross-section
  // in our definition inelastic and fission.
  
  theZ = Z;
  theA = A;
  theNames.SetMaxOffSet(5);
  G4NeutronHPDataUsed dataUsed;
  G4String rest = "/CrossSection/";
  G4String base = getenv("NeutronHPCrossSections");
  G4String base1 = base + "/Inelastic/";
  G4bool hasInelasticData;
  dataUsed = theNames.GetName(A, Z, base1, rest, hasInelasticData);
  G4String aName = dataUsed.GetName();
  G4NeutronHPVector inelasticData;
  G4double dummy;
  if(hasInelasticData)
  {
    ifstream aDataSet(aName, ios::in);
    aDataSet >> dummy >> dummy;
    inelasticData.Init(aDataSet, eV);
  }
  rest = "/CrossSection/";
  base1 = base + "/Fission/";
  G4bool hasFissionData = false;
  if(Z>=91)
  {
    dataUsed = theNames.GetName(A, Z, base1, rest, hasFissionData);
    aName = dataUsed.GetName();
  }
  G4NeutronHPVector fissionData;
  if(hasFissionData)
  {
    ifstream aDataSet(aName, ios::in);
    aDataSet >> dummy >> dummy;
    fissionData.Init(aDataSet, eV);
  }
  hasData = hasFissionData||hasInelasticData;
  if(hasData)
  {
    if(hasFissionData&&hasInelasticData)
    {
      theCrossSection.Merge(&fissionData, &inelasticData);
    }
    else if(hasFissionData)
    {
      theCrossSection = fissionData;
    }
    else if(hasInelasticData)
    {
      theCrossSection = inelasticData;
    }
    theCrossSection.Times(frac);
  }
  
  // now isotope-production cross-sections
  theNames.SetMaxOffSet(1);
  rest = "/CrossSection/";
  base1 = base + "/IsotopeProduction/";
  G4bool hasIsotopeProductionData;
  dataUsed = theNames.GetName(A, Z, base1, rest, hasIsotopeProductionData);
  aName = dataUsed.GetName();
  if(hasIsotopeProductionData)
  {
    ifstream aDataSet(aName, ios::in);
    aDataSet>>theNumberOfProducts;
    theProductionData = new G4NeutronIsoProdCrossSections * [theNumberOfProducts];
    for(G4int i=0; i<theNumberOfProducts; i++)
    {
      G4String aName;
      aDataSet >> aName;
      theProductionData[i] = new G4NeutronIsoProdCrossSections(aName);
      theProductionData[i]->Init(aDataSet);
    }
  }
  else
  {
    hasData = false;
  }
}

G4double G4NeutronIsoIsoCrossSections::
GetCrossSection(G4double anEnergy)
{
  G4double result;
  result = theCrossSection.GetY(anEnergy);
  return result;
}

G4String G4NeutronIsoIsoCrossSections::
GetProductIsotope(G4double anEnergy)
{
  G4String result;
  G4double * xSec = new G4double[theNumberOfProducts];
  G4double sum = 0;
  {
  for(G4int i=0; i<theNumberOfProducts; i++)
  {
    xSec[i] = theProductionData[i]->GetProductionCrossSection(anEnergy);
    sum += xSec[i];
  }
  }
  G4double random = G4UniformRand();
  G4double running = 0;
  G4int index;
  {
  for(G4int i=0; i<theNumberOfProducts; i++)
  {
    running += xSec[i];
    index = i;
    if(random<=running/sum) break;
  }
  }
  delete [] xSec;
  result = theProductionData[index]->GetProductIsotope();
  
  return result;
}
