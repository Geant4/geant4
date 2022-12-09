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
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
//	History: first implementation
//      June 2017 K.L. Genser, K. Lynch
//
// ---------------------------------------------------------------
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4MuonicAtomDecayPhysics.hh"

#include "G4GenericMuonicAtom.hh"
#include "G4MuonicAtomDecay.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"
#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4MuonicAtomDecayPhysics);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuonicAtomDecayPhysics::G4MuonicAtomDecayPhysics(G4int vL)
   :  G4VPhysicsConstructor("G4MuonicAtomDecay")
{
  G4PhysListUtil::InitialiseParameters();
  SetVerboseLevel(vL);
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
   G4cout << "G4MuonicAtomDecayPhysics() with verboseLevel "
          << verboseLevel << G4endl;
  }
#endif
}

G4MuonicAtomDecayPhysics::G4MuonicAtomDecayPhysics(const G4String&)
  : G4MuonicAtomDecayPhysics(1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MuonicAtomDecayPhysics::~G4MuonicAtomDecayPhysics()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuonicAtomDecayPhysics::ConstructParticle()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout << "G4MuonicAtomDecayPhysics::ConstructParticle()" << G4endl;
  }
#endif
  G4GenericMuonicAtom::GenericMuonicAtom();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MuonicAtomDecayPhysics::ConstructProcess()
{
  G4PhysicsListHelper* plh = G4PhysicsListHelper::GetPhysicsListHelper();
  G4int plhvl = plh->GetVerboseLevel();
#ifdef G4VERBOSE
  if (GetVerboseLevel()>0) {
    G4cout << "G4MuonicAtomDecayPhysics::ConstructProcess() verboseLevel "
           << GetVerboseLevel() << " to be set to " << plhvl << G4endl;
  }
#endif
  G4bool rc1 = plh->
    RegisterProcess(new G4MuonicAtomDecay(nullptr), // default G4HadronicInteraction* will be set
                    G4GenericMuonicAtom::GenericMuonicAtom());
  plh->SetVerboseLevel(plhvl);
  if (!(rc1)) {
    G4cout << " G4MuonicAtomDecayPhysics::ConstructProcess() : "
          << " RegisterProcess failed for G4GenericMuonicAtom "
          << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

