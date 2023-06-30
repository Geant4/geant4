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
// Standalone test of the G4 <-> FLUKA interface.
// One can simply set the hadronic interaction of interest,
// and gets in return the XS and final state.
//
// Note that this is just a test and NOT A PROPER G4 APPLICATION
// (no UI, no G4 event loop, NO CALL TO run/initialize etc).
//
// It relies on the fact that the physics processes are created 
// at physics list construction time anyway,
// and that the physics initialization happens 
// with runManager->SetUserInitialization(physicsList).
//
// Author: G.Hugo, 01 August 2022
//
// ***************************************************************************
#ifdef G4_USE_FLUKA


// G4
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"
//#include "G4ParticleDefinition.hh"
#include "G4HadProjectile.hh"
// G4
//#include "G4PionMinus.hh"
#include "G4Proton.hh"
// G4
#include "G4RunManagerFactory.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"

#include "fluka_interface.hh"
#include "FLUKAParticleTable.hh"
#include "G4_HP_CernFLUKAHadronInelastic_PhysicsList.hh"


G4int main() {

  G4cout << "Starting test_interface" << G4endl;

  // Construct a serial RUN MANAGER.
  std::unique_ptr<G4RunManager> runManager(G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly));
        
  // Create physics list with hadron inelastic processes from FLUKA.
  G4VModularPhysicsList* physicsList = new G4_HP_CernFLUKAHadronInelastic_PhysicsList();
  runManager->SetUserInitialization(physicsList);
  fluka_particle_table::initialize();
  //const auto& ionTable = G4ParticleTable::GetParticleTable()->GetIonTable();


  // TUNE HERE: PARTICLE
  //const auto particleKind = G4PionMinus::PionMinus();
  const auto particleKind = G4Proton::Proton();

  // TUNE HERE: KINETIC ENERGY
  const G4int numEvents = 1;
  const G4double minKineticEnergy = 1. * CLHEP::GeV;
  const G4double maxKineticEnergy = 1. * CLHEP::GeV;

  const G4double deltaKineticEnergy = (numEvents == 1 ? 0.
                                     : (maxKineticEnergy - minKineticEnergy) / (numEvents - 1));
  G4double kineticEnergy = minKineticEnergy;

  // TUNE HERE: DIRECTION
  const G4ThreeVector momentumDirection = G4ThreeVector(0., 0., 1.);
	
  // TUNE HERE: TARGET MATERIAL
  const G4int targetA = 28;
  const G4int targetZ = 14;
  const G4Nucleus targetNucleus = G4Nucleus(targetA, targetZ);
		
  G4HadFinalState finalState = G4HadFinalState();

  // PRINTOUT REQUESTED INTERACTION(S)
  G4cout << "particle = " << particleKind->GetParticleName() << G4endl;
  G4cout << "numEvents = " << numEvents << G4endl;
  G4cout << "minKineticEnergy = " << minKineticEnergy << G4endl;
  G4cout << "maxKineticEnergy = " << maxKineticEnergy << G4endl;
  G4cout << "deltaKineticEnergy = " << deltaKineticEnergy << G4endl;
  G4cout << "momentumDirection = " << momentumDirection << G4endl;
  G4cout << "targetA = " << targetA << G4endl;
  G4cout << "targetZ = " << targetZ << G4endl;
			   

  // Event loop		
  for (G4int eventIndex = 0; eventIndex < numEvents; ++eventIndex) {

    G4cout << "Start event index = " << eventIndex << G4endl;

    const G4DynamicParticle dynamicParticle = G4DynamicParticle(particleKind, 
                                                                momentumDirection,
                                                                kineticEnergy);

    // XS
    const G4double inelasticXS = fluka_interface::computeInelasticScatteringXS(&dynamicParticle, 
                                                                             targetZ, 
                                                                             targetA);

    G4cout << "fluka_interface::computeInelasticScatteringXS, inelasticXS = " << inelasticXS/barn << " [barn]." << G4endl;

    // FS
    finalState.Clear();
    finalState.SetStatusChange(stopAndKill);
    const G4HadProjectile projectile = G4HadProjectile(dynamicParticle);

    fluka_interface::setNuclearInelasticFinalState(&finalState,
                                                   projectile,
                                                   targetNucleus);

    // Update kinetic energy for the next event
    kineticEnergy += deltaKineticEnergy;
  }

}


#endif // G4_USE_FLUKA
