#include "G4AlphaPHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"

G4AlphaPHPBuilder::
G4AlphaPHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
}

G4AlphaPHPBuilder::
~G4AlphaPHPBuilder() 
{
  delete theParticlePHPModel;
}

void G4AlphaPHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4AlphaPHPBuilder::Build() not adding elastic" << G4endl;
}

void G4AlphaPHPBuilder::
Build(G4AlphaInelasticProcess * aP)
{
  G4cout << " G4AlphaPHPBuilder " << G4endl;
  G4ParticleHPInelasticData* theAlphaHPInelasticData=new G4ParticleHPInelasticData(G4Alpha::Alpha());
  theAlphaHPInelasticData->SetMinKinEnergy(theMin);
  theAlphaHPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theAlphaHPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4Alpha::Alpha(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);

}

