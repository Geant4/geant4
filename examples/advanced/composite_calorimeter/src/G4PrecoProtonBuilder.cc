#include "G4PrecoProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4PrecoProtonBuilder::
G4PrecoProtonBuilder() 
{
  theMin = 0;
  theMax = 150.*MeV;
  theModel = new G4PreCompoundModel(&theHandler);
}

G4PrecoProtonBuilder::
~G4PrecoProtonBuilder() {}

void G4PrecoProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4PrecoProtonBuilder::
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
