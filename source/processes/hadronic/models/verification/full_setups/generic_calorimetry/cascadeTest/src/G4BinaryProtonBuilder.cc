#include "G4BinaryProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BinaryProtonBuilder::
G4BinaryProtonBuilder() 
{
  theMin = 0;
  theModel = new G4BinaryCascade;
}

G4BinaryProtonBuilder::
~G4BinaryProtonBuilder() {}

void G4BinaryProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4BinaryProtonBuilder::
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
