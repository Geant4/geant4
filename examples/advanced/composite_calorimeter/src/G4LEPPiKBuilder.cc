#include "G4LEPPiKBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

G4LEPPiKBuilder::
G4LEPPiKBuilder()  
{
  theElasticModel = new G4LElastic();
  theMin = 0;
  theMax = 55*GeV;
}

G4LEPPiKBuilder::
~G4LEPPiKBuilder() {}

void G4LEPPiKBuilder::
Build(G4HadronElasticProcess & aP)
{
  aP.RegisterMe(theElasticModel);
}

void G4LEPPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  theLEPionPlusModel = new G4LEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(theMin);
  theLEPionPlusModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theLEPionPlusModel);
}

void G4LEPPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  theLEPionMinusModel = new G4LEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(theMin);
  theLEPionMinusModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theLEPionMinusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  theLEKaonPlusModel->SetMinEnergy(theMin);
  theLEKaonPlusModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theLEKaonPlusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  theLEKaonMinusModel->SetMaxEnergy(theMax);
  theLEKaonMinusModel->SetMinEnergy(theMin);
  aP.RegisterMe(theLEKaonMinusModel);
}

void G4LEPPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theLEKaonZeroLModel->SetMaxEnergy(theMax);
  theLEKaonZeroLModel->SetMinEnergy(theMin);
  aP.RegisterMe(theLEKaonZeroLModel);
}
 
void G4LEPPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theLEKaonZeroSModel->SetMaxEnergy(theMax);
  theLEKaonZeroSModel->SetMinEnergy(theMin);
  aP.RegisterMe(theLEKaonZeroSModel);
}

// 2002 by J.P. Wellisch
