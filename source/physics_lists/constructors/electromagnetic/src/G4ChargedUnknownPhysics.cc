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
// --------------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4ChargedUnknownPhysics.cc
//
//      Author:        A.Ribon
// 
//      Creation date: August 2024
//
//      Description:   This physics list constructor class is similar to
//                     G4UnknownDecayPhysics (which assigns the decay process
//                     (and transportation as well) to unknown particles),
//                     for charged unknown particles: it assigns the two EM
//                     processes of ionisation and multiple scattering
//                     (and transportation as well).
//
//      Modifications:
//
// --------------------------------------------------------------------------
//

#include "G4ChargedUnknownPhysics.hh"
#include "G4ChargedUnknownParticle.hh"
#include "G4UnknownParticle.hh"
#include "G4ProcessManager.hh"
#include "G4DynamicParticleMSC.hh"
#include "G4DynamicParticleIonisation.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY( G4ChargedUnknownPhysics );


G4ChargedUnknownPhysics::G4ChargedUnknownPhysics( G4int ver )
  : G4VPhysicsConstructor( "ChargedUnknownPhysics" ), verbose( ver ) {}


G4ChargedUnknownPhysics::G4ChargedUnknownPhysics( const G4String& name, G4int ver )
  : G4VPhysicsConstructor( name ), verbose( ver ) {}


void G4ChargedUnknownPhysics::ConstructParticle() {
  // Although this class deals only with 'chargedunknown' particles,
  // it is more consistent to define also 'unknown' particles.
  G4ChargedUnknownParticle::ChargedUnknownParticleDefinition();
  G4UnknownParticle::UnknownParticleDefinition();
}


void G4ChargedUnknownPhysics::ConstructProcess() {
  G4ProcessManager* pmanager = G4ChargedUnknownParticle::Definition()->GetProcessManager();
  pmanager->AddProcess( new G4DynamicParticleMSC, -1, 1, -1 );
  pmanager->AddProcess( new G4DynamicParticleIonisation, -1, 2, 2 );
  if ( verbose > 1 ) {
    G4cout << "### Added ionization and multiple scattering for  chargedunknown" << G4endl;
  }
}
