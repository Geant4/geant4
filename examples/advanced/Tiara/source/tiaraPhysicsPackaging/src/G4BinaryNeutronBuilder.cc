#include "G4BinaryNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BinaryNeutronBuilder::
G4BinaryNeutronBuilder() 
{
  theMin = 0;
  theMax = 3*GeV;
  theModel = new G4BinaryCascade;
}

G4BinaryNeutronBuilder::
~G4BinaryNeutronBuilder() {}

void G4BinaryNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4BinaryNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4BinaryNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4BinaryNeutronBuilder::
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
