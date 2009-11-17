#include "ML2PhysicsList.h"


CML2PhysicsList::CML2PhysicsList(void):  G4VUserPhysicsList() 
{
	this->SetDefaultCutValue(1.0*mm);
	this->SetCutsWithDefault();
}

CML2PhysicsList::~CML2PhysicsList(void)
{}
void CML2PhysicsList::ConstructParticle()
{
	G4Electron::ElectronDefinition();
	G4Gamma::GammaDefinition();
	G4Positron::PositronDefinition();
}
void CML2PhysicsList::ConstructProcess()
{
	AddTransportation();
	ConstructEM();
}
void CML2PhysicsList::ConstructEM()
{
 // Add standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {

      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);

    } else if (particleName == "e+") {

      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);

	}
	}
}
void CML2PhysicsList::ShowCutsValues()
{
	G4Region *reg;
	G4String regName;
	G4double cutValue;
	G4RegionStore *regStore=G4RegionStore::GetInstance();
	G4int nRegions=regStore->size();
	std::vector <G4Region*>::iterator reg_Iter;
	reg_Iter=regStore->begin();
	for (int i=0; i< nRegions;i++)
	{
		reg=reg_Iter[i];
		regName = reg->GetName();
		cutValue=reg->GetProductionCuts()->GetProductionCut(0);
		std::cout << regName<<" cut Value: " << cutValue/mm<<" [mm]"<< G4endl;
	}
}
void CML2PhysicsList::SetCuts()
{
	this->ShowCutsValues();
}

