#include "G4QGSPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4QGSPNeutronBuilder::
G4QGSPNeutronBuilder() 
{
  theMin = 15*GeV;
  theModel = new G4TheoFSGenerator;
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
}

G4QGSPNeutronBuilder::
~G4QGSPNeutronBuilder() 
{
  delete theStringDecay;
}

void G4QGSPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4QGSPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
