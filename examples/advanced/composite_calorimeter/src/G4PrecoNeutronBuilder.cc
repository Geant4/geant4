#include "G4PrecoNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4PrecoNeutronBuilder::
G4PrecoNeutronBuilder() 
{
  theMin = 0;
  theMax = 150.*MeV;
  theModel = new G4PreCompoundModel(&theHandler);
}

G4PrecoNeutronBuilder::
~G4PrecoNeutronBuilder() {}

void G4PrecoNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4PrecoNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
