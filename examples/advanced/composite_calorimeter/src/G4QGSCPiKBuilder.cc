#include "G4QGSCPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4QGSCPiKBuilder::
G4QGSCPiKBuilder() 
{
  theMin = 15*GeV;
  theModel = new G4TheoFSGenerator;
  theCascade = new G4StringChipsParticleLevelInterface;
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
}

G4QGSCPiKBuilder::
~G4QGSCPiKBuilder() 
{
  delete theStringDecay;
}

void G4QGSCPiKBuilder::
Build(G4HadronElasticProcess & aP) {}

void G4QGSCPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

void G4QGSCPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

void G4QGSCPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

void G4QGSCPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

void G4QGSCPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

void G4QGSCPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
}

// 2002 by J.P. Wellisch
