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
// Base class to define a common interface for post-cascade processing.
//
// 20130621  Add implementation of ::collide() to throw exception

#include "G4VCascadeDeexcitation.hh"
#include "G4HadronicException.hh"


// Standard Collider interface must not be used, only G4Fragment

void G4VCascadeDeexcitation::collide(G4InuclParticle* /*bullet*/,
				     G4InuclParticle* /*target*/,
				     G4CollisionOutput& /*globalOutput*/) {
  if (verboseLevel) {
    G4cout << " >>> G4VCascadeDeexcitation[" << theName << "]::collide "
	   << " *** SHOULD NOT BE CALLED ***" << G4endl;
  }

  throw G4HadronicException(__FILE__, __LINE__, 
    "G4VCascadeDeexcitation::collide() invalid, must use ::deExcite(G4Fagment*)");
}
