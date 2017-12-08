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
// ClassName:   G4StepLimiterPhysics
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

#include "G4StepLimiterPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4StepLimiterPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterPhysics::G4StepLimiterPhysics(const G4String& name)
  :  G4VPhysicsConstructor(name),fApplyToAll(false)
{
  SetPhysicsType(bUnknown);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterPhysics::~G4StepLimiterPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4StepLimiterPhysics::ConstructParticle()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4StepLimiterPhysics::ConstructProcess()
{
  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();

  G4StepLimiter* stepLimiter = new G4StepLimiter();
  G4UserSpecialCuts* userSpecialCuts = new G4UserSpecialCuts();
  while ((*myParticleIterator)()) {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4double charge = particle->GetPDGCharge();

    if(!particle->IsShortLived()) {
      if (charge != 0.0 || fApplyToAll) {
	// All charged particles should have a step limiter
	// to make sure that the steps do not get too long.
	pmanager->AddDiscreteProcess(stepLimiter);
	pmanager->AddDiscreteProcess(userSpecialCuts);
      } else {
	// Energy cuts for all other neutral particles
	pmanager->AddDiscreteProcess(userSpecialCuts);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
