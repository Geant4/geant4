//    **********************************
//    *                                *
//    *     BrachyPhysicsList.cc       *
//    *                                *
//    **********************************

#include "BrachyPhysicsList.hh"

#include "G4EnergyLossTables.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MultipleScattering.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

//....

BrachyPhysicsList::BrachyPhysicsList()
{
 m_pPositron = NULL;
 m_pGamma = NULL;
 m_pElectron = NULL;
}

//....

BrachyPhysicsList::~BrachyPhysicsList()
{
}

//....

void BrachyPhysicsList::ConstructParticle()
{
 m_pGamma = G4Gamma::GammaDefinition();
 m_pElectron = G4Electron::ElectronDefinition();
 m_pPositron = G4Positron::PositronDefinition();
 
}

//....

void BrachyPhysicsList::ConstructProcess()
{
 AddTransportation();
 ConstructEM();
}

//....

void BrachyPhysicsList::ConstructEM()
{
 G4ProcessManager *pProcessManager = m_pElectron->GetProcessManager();
 if(pProcessManager)
	{
	pProcessManager->AddProcess(new G4MultipleScattering(),-1, 1,1);
      pProcessManager->AddProcess(new G4eIonisation(),       -1, 2,2);
      pProcessManager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
	}
 pProcessManager = m_pGamma->GetProcessManager();
 if(pProcessManager)
	{
	pProcessManager->AddDiscreteProcess(new G4PhotoElectricEffect());
	pProcessManager->AddDiscreteProcess(new G4ComptonScattering());
	pProcessManager->AddDiscreteProcess(new G4GammaConversion());
	}
 pProcessManager = m_pPositron->GetProcessManager();
 if(pProcessManager)
	{
      pProcessManager->AddProcess(new G4MultipleScattering(),-1, 1,1);
      pProcessManager->AddProcess(new G4eIonisation(),       -1, 2,2);
      pProcessManager->AddProcess(new G4eBremsstrahlung(),   -1,-1,3);
      pProcessManager->AddProcess(new G4eplusAnnihilation(),  0,-1,4);      
      }
}

//....

void BrachyPhysicsList::SetCuts()
{
 // uppress error messages even in case e/gamma/proton do not exist            
 G4int temp = GetVerboseLevel();
 SetVerboseLevel(0);                                                           
 SetCutsWithDefault();
 // Retrieve verbose level
 SetVerboseLevel(temp);  
}

