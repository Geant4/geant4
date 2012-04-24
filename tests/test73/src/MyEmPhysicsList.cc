
// $Id: MyEmPhysicsList.cc,v 1.3 2009-11-14 18:04:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "MyEmPhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4CoulombScattering.hh"
#include "G4WentzelVIModel.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyEmPhysicsList::MyEmPhysicsList(const G4String& name)
:  G4VPhysicsConstructor(name), verbose(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyEmPhysicsList::~MyEmPhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyEmPhysicsList::ConstructProcess()
{
	// Add standard EM Processes

	theParticleIterator->reset();
	while( (*theParticleIterator)() )
	{
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		if (particleName == "gamma")
		{
			// gamma         
			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
			pmanager->AddDiscreteProcess(new G4ComptonScattering);
			pmanager->AddDiscreteProcess(new G4GammaConversion);

		} 
		
		else if (particleName == "e-")
		{
			//electron
			G4eMultipleScattering* msc = new G4eMultipleScattering();
			//msc->AddEmModel(0, new G4UrbanMscModel93());
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,                   -1, 1, 1);
			G4eIonisation* eIoni = new G4eIonisation();
			eIoni->SetStepFunction(0.2, 100*um);      
			pmanager->AddProcess(eIoni,                 -1, 2, 2);
			pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);

		} 
		
		else if (particleName == "e+")
		{
			//positron
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
			pmanager->AddProcess(new G4eIonisation,         -1, 2,2);
			pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3,3);
			pmanager->AddProcess(new G4eplusAnnihilation,    0,-1,4);

		}
		
		else if( particleName == "mu+" || 
				particleName == "mu-"    ) 
		{
			//muon  
			G4MuMultipleScattering* msc = new G4MuMultipleScattering();
			msc->AddEmModel(0, new G4WentzelVIModel());
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,                       -1, 1, 1);
			G4MuIonisation* muIoni = new G4MuIonisation();
			muIoni->SetStepFunction(0.2, 50*um);          
			pmanager->AddProcess(muIoni,                    -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
			pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);
			//AddStepMax(particle, pmanager);
			//pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5);
			//pmanager->AddDiscreteProcess(new G4CoulombScattering());
		}
		
		else if( particleName == "alpha" || particleName == "GenericIon" ) 
		{ 
			pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
			pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
		}

		else if (particleName == "pi+" ||
				particleName == "pi-" ||
				particleName == "kaon+" ||
				particleName == "kaon-" ||
				particleName == "proton" ) 
		{

			//pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			G4hMultipleScattering* msc = new G4hMultipleScattering();
			msc->AddEmModel(0, new G4WentzelVIModel());
			msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,                   -1, 1, 1);

			G4hIonisation* hIoni = new G4hIonisation();
			hIoni->SetStepFunction(0.2, 50*um);
			pmanager->AddProcess(hIoni,                     -1, 2, 2);
			pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
			pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
			//pmanager->AddDiscreteProcess(new G4CoulombScattering());
		}

		else if ((!particle->IsShortLived()) &&
				(particle->GetPDGCharge() != 0.0) && 
				(particle->GetParticleName() != "chargedgeantino")) 
		{
			//all others charged particles except geantino
			pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
			pmanager->AddProcess(new G4hIonisation,        -1,2,2);
		}

		G4EmProcessOptions opt;
		opt.SetVerbose(verbose);
		//opt.SetPolarAngleLimit(0.2);
		opt.SetPolarAngleLimit(CLHEP::pi);
		opt.SetApplyCuts(true);

	}
}
// To limit step size in logical volumes set in Detector Geometry class
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void MyEmPhysicsList::AddStepMax(G4ParticleDefinition* particle,
		G4ProcessManager* pmanager)
{
	// Step limitation seen as a process
	G4StepLimiter* stepLimiter = new G4StepLimiter();
	////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();

		if (particle->GetPDGCharge() != 0.0)
		{
			pmanager ->AddDiscreteProcess(stepLimiter);
			////pmanager ->AddDiscreteProcess(userCuts);
		}
	
}

	//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

