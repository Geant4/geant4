#include "G4LEADProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEADProtonBuilder::
G4LEADProtonBuilder() 
{
  theMin = 0;
  theModel = new G4Mars5GeV;
}

G4LEADProtonBuilder::
~G4LEADProtonBuilder() {}

void G4LEADProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4LEADProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * thePStore;
  thePStore = aP.GetCrossSectionDataStore();
  thePStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
