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
//
// $Id: G4EmHadronBuilder.cc,v 1.2 2006/06/29 17:27:21 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

