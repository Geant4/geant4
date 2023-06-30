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
#ifdef G4_USE_FLUKA


#include "FLUKANuclearInelasticModel.hh"

// G4
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"

#include "fluka_interface.hh"


// ***************************************************************************
// FLUKA hadron inelastic physics final state.
// ***************************************************************************
FLUKANuclearInelasticModel::FLUKANuclearInelasticModel()
  : G4HadronicInteraction("FLUKANuclearInelasticModel"),
    finalState_(std::make_unique<G4HadFinalState>())
{}


// ***************************************************************************
// FLUKA hadron inelastic physics: returns final state from FLUKA.
// ***************************************************************************
G4HadFinalState* FLUKANuclearInelasticModel::ApplyYourself(const G4HadProjectile& projectile, 
							   G4Nucleus& targetNucleus) {

  // Clean-up final state.
  finalState_->Clear();
  finalState_->SetStatusChange(stopAndKill);
	
  // GET FINAL STATE FROM FLUKA INTERFACE
  fluka_interface::setNuclearInelasticFinalState(finalState_.get(),
                                                 projectile,
                                                 targetNucleus);

  return finalState_.get();
}


#endif //G4_USE_FLUKA
