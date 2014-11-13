#include "G4ProtonPHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"

G4ProtonPHPBuilder::
G4ProtonPHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
}

G4ProtonPHPBuilder::
~G4ProtonPHPBuilder() 
{
  delete theParticlePHPModel;
}

void G4ProtonPHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4ProtonPHPBuilder::Build() not adding elastic" << G4endl;
}

void G4ProtonPHPBuilder::
Build(G4ProtonInelasticProcess * aP)
{
  G4ParticleHPInelasticData* theProtonHPInelasticData=new G4ParticleHPInelasticData(G4Proton::Proton());
  theProtonHPInelasticData->SetMinKinEnergy(theMin);
  theProtonHPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theProtonHPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4Proton::Proton(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);
  //  theParticleModel->DumpInfo();

}

