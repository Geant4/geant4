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


#include "G4_HP_CernFLUKAHadronInelastic_PhysicsList.hh"

// G4
//#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"         // LIV
#include "G4EmExtraPhysics.hh"
// G4
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"    // HP
// G4
//#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"     // HP
// G4
//#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "FLUKAHadronInelasticPhysicsConstructor.hh"
// G4
#include "G4IonPhysics.hh"
// G4
#include "G4StoppingPhysics.hh"
//#include "G4NeutronTrackingCut.hh"        // when NOT HP
// G4
#include "CLHEP/Units/SystemOfUnits.h"

#include "fluka_interface.hh"


// ***************************************************************************
// A physics list based on FTFP_BERT_HP LIV, 
// but with full hadron inelastic physics replaced by the one from FLUKA.CERN.
// ***************************************************************************
G4_HP_CernFLUKAHadronInelastic_PhysicsList::G4_HP_CernFLUKAHadronInelastic_PhysicsList(G4int verbose) {

  if (verbose > 0) {
    G4cout << "<<< Geant4 Physics List simulation engine: G4_HP_CernFLUKAHadronInelastic_PhysicsList" 
           << G4endl << G4endl;
  }
  SetVerboseLevel(verbose);

  // IMPORTANT: The default production cut is set here.
  defaultCutValue = 0.7*CLHEP::mm;
  

  // EM PHYSICS
  //RegisterPhysics( new G4EmStandardPhysics(ver) );
  RegisterPhysics( new G4EmLivermorePhysics( verbose ) );    // LIV

  // SYNCHROTON RADIATION & GN PHYSICS
  RegisterPhysics( new G4EmExtraPhysics(verbose) );

  // DECAYS 
  RegisterPhysics( new G4DecayPhysics(verbose) );
  RegisterPhysics( new G4RadioactiveDecayPhysics(verbose) ); // HP

  // HADRON ELASTIC SCATTERING
  //RegisterPhysics( new G4HadronElasticPhysics(verbose) );
  RegisterPhysics( new G4HadronElasticPhysicsHP(verbose) );   // HP

  // HADRON INELASTIC PHYSICS
  //RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(verbose) );
  RegisterPhysics( new FLUKAHadronInelasticPhysicsConstructor( verbose ) );

  // ION PHYSICS
  RegisterPhysics( new G4IonPhysics(verbose) );

  // STOPPING PHYSICS
  RegisterPhysics( new G4StoppingPhysics(verbose) );
  // NEUTRON TRACKING CUT
  // NB: Not in FTFP_BERT_HP!
  //RegisterPhysics( new G4NeutronTrackingCut( verbose ) ); // when NOT HP

  
  // IMPORTANT: Initialize the FLUKA interface here.
  // Both activation switches should be set to TRUE to provide the most comprehensive results.
  // NB: COMPARISON WITH G4 DOES NOT SEEM MEANINGFUL 
  // WHEN COALESCENCE IS ACTIVATED IN BOTH FLUKA AND G4.
  // Freedom to choose & see the effect of these switches is hence provided here.
  const G4bool activateCoalescence = true;
  const G4bool activateHeavyFragmentsEvaporation = true;
  fluka_interface::initialize(activateCoalescence, 
                              activateHeavyFragmentsEvaporation);
}


// ***************************************************************************
// IMPORTANT: Set production cuts here: add a 0. cut for proton production.
// ***************************************************************************
void G4_HP_CernFLUKAHadronInelastic_PhysicsList::SetCuts() {

	if (verboseLevel > 1) { G4cout << "G4_HP_CernFLUKAHadronInelastic_PhysicsList::SetCuts:"; }

	// G4VUserPhysicsList::SetCutsWithDefault sets 
	// the default cut value for all particle types.
	SetCutsWithDefault();

	// Set proton cut value to 0, for producing low energy recoil nucleus.
	SetCutValue(0., "proton");
	G4cout << "Proton production cut: " << GetCutValue("proton") << G4endl;
}


#endif // G4_USE_FLUKA

