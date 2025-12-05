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
//----------------------------------------------------------------------------//
// This physics list *must* be set with a *reference* physics list.
//   These provide a full set of models (both electromagnetic and hadronic).
//
//   The reference physics list must be set by issuing its command in a macro
//   file (see messenger class). No other action is required.
//   Examples of physics list names: QGSP_BIC_HP_EMZ or QGSP_BERT_HP
//----------------------------------------------------------------------------//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  SetDefaultCutValue(1.0*mm);

  // This is to force setting the physics list with macro command
  fPhysListIsSet = false;

  fVerbose = 1;
  fMessenger = new PhysicsListMessenger(this);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Construct Particles" << G4endl;
  }
  // This method is invoked when the Geant4 application starts
  // (do not mix with run initialization).

  // (Taken from G4DecayPhysics)
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  if(fVerbose > 0) {
    G4cout << "### PhysicsList Construct Processes" << G4endl;
  }

  if (fPhysListIsSet)
    G4VModularPhysicsList::ConstructProcess();
  else
    G4Exception("PhysicsList::ConstructProcess()", "PhysList001",
		FatalException, "No PHYSICS LIST has been set!");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetPhysicsList(const G4String& name)
{
  if(fVerbose > 0)
    G4cout << "### PhysicsList set physics list <" << name
           << "> " << G4endl;

  if (!fPhysListIsSet) {
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = factory.GetReferencePhysList(name);

    size_t ii = 0;
    const G4VPhysicsConstructor* elem = phys->GetPhysics(ii);
    G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
    while (elem) {
      RegisterPhysics(tmp);
      G4cout << "PhysicsList Type: " << elem->GetPhysicsType() << G4endl;
      G4cout << "PhysicsList Name: " << elem->GetPhysicsName() << G4endl;
      elem = phys->GetPhysics(++ii);
      tmp = const_cast<G4VPhysicsConstructor*> (elem);
    }

    G4cout << name << " reference physics List has been ACTIVATED."
	   << G4endl;

    // Update the flag, the physics list is set
    fPhysListIsSet = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
