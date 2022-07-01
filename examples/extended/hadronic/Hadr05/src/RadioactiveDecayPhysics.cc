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

#include "RadioactiveDecayPhysics.hh"

#include "G4Radioactivation.hh"
#include "G4GenericIon.hh"
#include "G4PhysicsListHelper.hh"
#include "G4EmParameters.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclideTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RadioactiveDecayPhysics::RadioactiveDecayPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{
  // hadronic physics extra configuration
  //
  G4DeexPrecoParameters* deex = 
     G4NuclearLevelData::GetInstance()->GetParameters();     
  deex->SetStoreICLevelData(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                       /std::log(2.));
  deex->SetIsomerProduction(true);
  deex->SetCorrelatedGamma(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RadioactiveDecayPhysics::~RadioactiveDecayPhysics()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RadioactiveDecayPhysics::ConstructParticle()
{
  G4GenericIon::GenericIon();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RadioactiveDecayPhysics::ConstructProcess()
{
  G4Radioactivation* radioactiveDecay = new G4Radioactivation();

  G4bool ARMflag = false;
  radioactiveDecay->SetARM(ARMflag);        //Atomic Rearangement

  // EM physics extra configuration
  // this physics constructor should be defined after EM constructor
  G4EmParameters::Instance()->SetFluo(ARMflag);  
  G4EmParameters::Instance()->SetAugerCascade(ARMflag);
  G4EmParameters::Instance()->SetDeexcitationIgnoreCut(ARMflag);

  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();

  // EM physics constructors are not used
  if( ad == nullptr ) {
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    man->ResetParameters();
  }

  G4PhysicsListHelper::GetPhysicsListHelper()->
    RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

