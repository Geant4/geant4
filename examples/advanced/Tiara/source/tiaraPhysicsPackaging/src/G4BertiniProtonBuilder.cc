#include "G4BertiniProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BertiniProtonBuilder::
G4BertiniProtonBuilder() 
{
  theMin = 0;
  theModel = new G4CascadeInterface;
}

G4BertiniProtonBuilder::
~G4BertiniProtonBuilder() {}

void G4BertiniProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4BertiniProtonBuilder::
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
