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
// $Id: G4EmHadronBuilder.cc,v 1.1 2004/12/13 16:38:54 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-05 $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmHadronBuilder
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4EmHadronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4MultipleScattering.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHadronBuilder::G4EmHadronBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmHadronBuilder::~G4EmHadronBuilder()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHadronBuilder::ConstructParticle()
{
  G4Electron::Electron();
  G4Proton::Proton();
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmHadronBuilder::ConstructProcess()
{
  // Add standard EM Processes
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();

    if(particle->GetPDGMass() > 110.*MeV) {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "GenericIon" || particleName == "alpha" || particleName == "He3") {

        pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
        pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);

      } else if ((!particle->IsShortLived()) &&
	         (particle->GetPDGCharge() != 0.0) &&
	         (particle->GetParticleName() != "chargedgeantino")) {

        pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
        pmanager->AddProcess(new G4hIonisation,       -1,2,2);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

