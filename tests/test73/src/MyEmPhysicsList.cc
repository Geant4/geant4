//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "MyEmPhysicsList.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4UrbanMscModel93.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4MuMultipleScattering.hh"



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

	aParticleIterator->reset();
	while( (*aParticleIterator)() )
	{
		G4ParticleDefinition* particle = aParticleIterator->value();
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
            //msc->SetSkin(10);
			msc->AddEmModel(0, new G4UrbanMscModel93());
			msc->SetStepLimitType(fUseDistanceToBoundary);
            msc->SetSkin(10);
            G4cout<<"Using MSC model: G4UrbanMscModel93 with UseDistanceToBoundary steplimit type and Skin=10"<<G4endl;
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
			//msc->SetStepLimitType(fUseDistanceToBoundary);
			pmanager->AddProcess(msc,                       -1, 1, 1);
			G4MuIonisation* muIoni = new G4MuIonisation();
			muIoni->SetStepFunction(0.2, 50*um);          
			pmanager->AddProcess(muIoni,                    -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
			pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);
			//AddStepMax(particle, pmanager);
			//pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5);
			pmanager->AddDiscreteProcess(new G4CoulombScattering());
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
			G4MuMultipleScattering* msc = new G4MuMultipleScattering();
			msc->AddEmModel(0, new G4WentzelVIModel());
			//msc->SetStepLimitType(fUseDistanceToBoundary);
			//msc->SetSkin(10);

			pmanager->AddProcess(msc,                   -1, 1, 1);

			G4hIonisation* hIoni = new G4hIonisation();
			//hIoni->SetStepFunction(0.2, 50*um);
			pmanager->AddProcess(hIoni,                     -1, 2, 2);
			pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
			pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
			// this line is uncommented only because hadron 
			// elastic is not used in this Physics List 
			pmanager->AddDiscreteProcess(new G4CoulombScattering());
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

