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
// ------------------------------------------------------------
//      Bertini Cascade unboundPN class implementation file
//
//      History: first implementation, inspired by G4Proton
//      17 Nov 2009:  Michael Kelsey
//	06 Apr 2010:  Do G4Ions initialization in ctor.
//	13 Apr 2010:  Per Kurashige, inherit from G4VShortLivedParticle.
//	06 May 2010:  Remove created particle from master table.
//      25 May 2012:  Add flags to suppress particle-table error message.
//	01 May 2013:  Remove G4ThreadLocal from static pointer.
// ----------------------------------------------------------------

#include "G4UnboundPN.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                          UNBOUNDPN                             ###
// ######################################################################
G4UnboundPN* G4UnboundPN::theInstance = 0;

//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation
//          2*Isospin       2*Isospin3       G-parity
//               type    lepton number  baryon number   PDG encoding
//             stable         lifetime    decay table
G4UnboundPN::G4UnboundPN()
  : G4VShortLivedParticle("unboundPN",
			  (proton_mass_c2+neutron_mass_c2), 0.0*MeV, +1.*eplus, 
			  2,       +1,       0,          
			  2,        0,       0,             
			  "nucleus",        0,      +2, 0, /* ? 100010020 */
			  true,       0.,    NULL) {}

G4UnboundPN* G4UnboundPN::Definition() {
  if (0 == theInstance) {
    theInstance = new G4UnboundPN;	// There can be only one

    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4bool tableReady = pTable->GetReadiness();
    pTable->SetReadiness(false);        // Suppress error message
    pTable->Remove(theInstance);        // Make invisible to GEANT4
    pTable->SetReadiness(tableReady);   // Set back 'ready to use' flag
  }

  return theInstance;
}

// Simple call-throughs
G4UnboundPN* G4UnboundPN::UnboundPNDefinition() { return Definition(); }
G4UnboundPN* G4UnboundPN::UnboundPN()           { return Definition(); }
