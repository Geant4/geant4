//    **********************************
//    *                                *
//    *     BrachyPhysicsList.cc       *
//    *                                *
//    **********************************

#include "BrachyPhysicsList.hh"

// Comment this if you do not want to use the 
// LowEnergy models but prefer the standard
#define BRACHY_OPT_USELOWENERGY

#include "G4EnergyLossTables.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

#ifdef BRACHY_OPT_USELOWENERGY
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyRayleigh.hh"
#endif

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
 G4ProcessManager *pProcessManager;

 #ifdef BRACHY_OPT_USELOWENERGY	// LowEnergy models for e+, e- and gammas
 pProcessManager = m_pElectron->GetProcessManager();
 if(pProcessManager)
	{
	pProcessManager->AddProcess(new G4MultipleScattering(),-1, 1,1);
	pProcessManager->AddProcess(new G4LowEnergyIonisation(),       -1, 2,2);
	pProcessManager->AddProcess(new G4LowEnergyBremsstrahlung(),   -1,-1,3);
	}
 pProcessManager = m_pGamma->GetProcessManager();
 if(pProcessManager)
	{
	pProcessManager->AddDiscreteProcess(new G4LowEnergyPhotoElectric());
	pProcessManager->AddDiscreteProcess(new G4LowEnergyCompton());
	pProcessManager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
	pProcessManager->AddDiscreteProcess(new G4LowEnergyRayleigh());
	}
 pProcessManager = m_pPositron->GetProcessManager();
 if(pProcessManager)
	{
	pProcessManager->AddProcess(new G4MultipleScattering(),     -1, 1,1);
	pProcessManager->AddProcess(new G4LowEnergyIonisation(),    -1, 2,2);
	pProcessManager->AddProcess(new G4LowEnergyBremsstrahlung(),-1,-1,3);
	pProcessManager->AddProcess(new G4eplusAnnihilation(),       0,-1,4);      
	}
 #else				// Standard models for e+, e- and gammas
 pProcessManager = m_pElectron->GetProcessManager();
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
 #endif // BRACHY_OPT_USELOWENERGY
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

