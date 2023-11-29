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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RadioactiveDecayPhysics.hh"

#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"
#include "G4EmParameters.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"
#include "G4PhysListUtil.hh"
#include "G4Triton.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4RadioactiveDecayPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics::G4RadioactiveDecayPhysics(G4int ver)
:  G4VPhysicsConstructor("G4RadioactiveDecay")
{
  G4PhysListUtil::InitialiseParameters();
  SetVerboseLevel(ver);

  // hadronic physics extra configuration
  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetStoreICLevelData(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                       /std::log(2.));
  deex->SetIsomerProduction(true);
}

G4RadioactiveDecayPhysics::G4RadioactiveDecayPhysics(const G4String&, G4int ver)
  : G4RadioactiveDecayPhysics(ver)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics::~G4RadioactiveDecayPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics::ConstructProcess()
{
  // EM physics extra configuration
  // this physics constructor should be defined after EM constructor
  G4EmParameters::Instance()->SetAuger(true);
  G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);

  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();

  // EM physics constructors are not used
  if( ad == nullptr ) {
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    man->ResetParameters();
  }

  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(new G4RadioactiveDecay(), G4GenericIon::GenericIon());

  // Triton (which is not a generic ion) is the only light ion that decays.
  // Note that the anti_triton does not have beta decay, because RadioactiveDecay,
  // in its current implementation, does not handle any kind of anti-ions: 
  // in practice, this is an acceptable approximation because of its relatively
  // long lifetime and the fact that annihilation and nuclear capture
  // are more likely to happen before decay.
  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(new G4RadioactiveDecay(), G4Triton::Triton());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

