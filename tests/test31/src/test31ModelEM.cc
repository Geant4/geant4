//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// 

#include "test31ModelEM.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScatteringSTD.hh"

#include "G4eIonisationSTD.hh"
#include "G4eBremsstrahlungSTD.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisationSTD.hh"
#include "G4MuBremsstrahlungSTD.hh"
#include "G4MuPairProductionSTD.hh"

#include "G4hIonisationSTD.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31ModelEM::test31ModelEM(const G4String& name)
  :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31ModelEM::~test31ModelEM()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31ModelEM::ConstructParticle()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31ModelEM::ConstructProcess()
{
  // Add EM processes realised on base of prototype of model approach design
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      pmanager->AddProcess(new G4eIonisationSTD,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlungSTD,    -1,-1,3);

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      pmanager->AddProcess(new G4eIonisationSTD,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlungSTD,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,      0,-1,4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      //muon
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisationSTD,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlungSTD,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProductionSTD,  -1,-1,4);

    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      G4VEnergyLossSTD* ion = new G4ionIonisation();
      pmanager->AddProcess(ion,      -1,2,2);

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      G4VEnergyLossSTD* hIon = new G4hIonisationSTD();
      pmanager->AddProcess(hIon,   -1,2,2);
      //  hIon->SetLossFluctuations(false);
      hIon->SetDEDXBinning(240);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

