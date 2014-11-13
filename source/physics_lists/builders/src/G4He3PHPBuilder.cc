#include "G4He3PHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"

G4He3PHPBuilder::
G4He3PHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
}

G4He3PHPBuilder::
~G4He3PHPBuilder() 
{
  delete theParticlePHPModel;
}

void G4He3PHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4He3PHPBuilder::Build() not adding elastic" << G4endl;
}

void G4He3PHPBuilder::
Build(G4He3InelasticProcess * aP)
{
  G4ParticleHPInelasticData* theHe3HPInelasticData=new G4ParticleHPInelasticData(G4He3::He3());
  theHe3HPInelasticData->SetMinKinEnergy(theMin);
  theHe3HPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theHe3HPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4He3::He3(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);

}

