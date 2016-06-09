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
//---------------------------------------------------------------------------
//
// ClassName:   G4StepLimiterBuilder
//
// Author:      V.Ivanchenko 24.11.2004
//
// Modified:
//
//----------------------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4StepLimiterBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4StepLimiterBuilder);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterBuilder::G4StepLimiterBuilder(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  fStepLimiter = new G4StepLimiter();
  fUserSpecialCuts = new G4UserSpecialCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterBuilder::~G4StepLimiterBuilder()
{
  delete fStepLimiter;
  delete fUserSpecialCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4StepLimiterBuilder::ConstructParticle()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4StepLimiterBuilder::ConstructProcess()
{
  theParticleIterator->reset();

  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4double charge = particle->GetPDGCharge();

    if(!particle->IsShortLived()) {
      if (charge != 0.0) {
	// All charged particles should have a step limiter
	// to make sure that the steps do not get too long.
	pmanager->AddDiscreteProcess(fStepLimiter);
	pmanager->AddDiscreteProcess(fUserSpecialCuts);
      } else {
	// Energy cuts for all other neutral particles
	pmanager->AddDiscreteProcess(fUserSpecialCuts);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
