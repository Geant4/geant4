#include "G4CascadeProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4CascadeProtonBuilder::
G4CascadeProtonBuilder() 
{
  theMin = 0;
  theMax = 150.*MeV;
  theModel = new G4BinaryCascade();
}

G4CascadeProtonBuilder::
~G4CascadeProtonBuilder() {}

void G4CascadeProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4CascadeProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * thePStore;
  thePStore = aP.GetCrossSectionDataStore();
  thePStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
