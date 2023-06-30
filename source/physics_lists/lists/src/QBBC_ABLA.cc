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
//---------------------------------------------------------------------------
//
// ClassName:  QBBC_ABLA
//
// Author: Alberto Ribon (CERN), April 2023
//
// The new, experimental physics list QBBC_ABLA is similar to the reference
// physics list QBBC, except that for hadron inelastic the physics constructor
// G4HadronInelasticQBBC_ABLA is used (instead of G4HadronInelasticQBBC):
// in practice, QBBC_ABLA behaves as QBBC, with the only difference that for
// the final-state of nuclear inelastic interactions of charged pions and
// nucleons projectiles, the ABLA model (instead of the usual
// Precompound/de-excitation) is utilized for nuclear de-excitation.
//
// Modified:
//
//----------------------------------------------------------------------------
//

#include "QBBC_ABLA.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronInelasticQBBC_ABLA.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"


QBBC_ABLA::QBBC_ABLA( G4int ver, const G4String& ) {
  if ( ver > 0 ) G4cout << "<<< Experimental Reference Physics List QBBC_ABLA " << G4endl;
  defaultCutValue = 0.7*CLHEP::mm;
  SetVerboseLevel( ver );
  RegisterPhysics( new G4EmStandardPhysics(ver) );
  RegisterPhysics( new G4EmExtraPhysics(ver) );
  RegisterPhysics( new G4DecayPhysics(ver) );
  RegisterPhysics( new G4HadronElasticPhysicsXS(ver) );
  RegisterPhysics( new G4StoppingPhysics(ver) );
  RegisterPhysics( new G4IonPhysicsXS(ver) );
  RegisterPhysics( new G4IonElasticPhysics(ver) );
  RegisterPhysics( new G4HadronInelasticQBBC_ABLA(ver) );
  RegisterPhysics( new G4NeutronTrackingCut(ver) );
}		 
