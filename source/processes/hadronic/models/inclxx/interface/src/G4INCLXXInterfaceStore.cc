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
// INCL++ intra-nuclear cascade model
// Pekka Kaitaniemi, CEA and Helsinki Institute of Physics
// Davide Mancusi, CEA
// Alain Boudard, CEA
// Sylvie Leray, CEA
// Joseph Cugnon, University of Liege
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

/** \file G4INCLXXInterfaceStore.cc
 * \brief The G4INCLXXInterfaceStore class implementation
 *
 * \date 24 May 2012
 * \author Davide Mancusi
 */

#include "G4INCLXXInterfaceMessenger.hh"

G4INCLXXInterfaceStore *G4INCLXXInterfaceStore::theInstance = NULL;

G4INCLXXInterfaceStore::G4INCLXXInterfaceStore() :
  dumpInput(false),
  accurateProjectile(true),
  theMaxClusterMassDefault(8),
  theMaxClusterMass(theMaxClusterMassDefault),
  theMaxProjMassINCL(18),
  theINCLModel(NULL),
  nWarnings(0),
  maxWarnings(50)
{
  theINCLXXInterfaceMessenger = new G4INCLXXInterfaceMessenger(this);
}

G4INCLXXInterfaceStore::~G4INCLXXInterfaceStore() {
  delete theINCLXXInterfaceMessenger;
  delete theINCLModel;
}

void G4INCLXXInterfaceStore::EmitWarning(const G4String &message) {
  if(++nWarnings<=maxWarnings) {
    G4cout << "[INCL++] Warning: " << message << G4endl;
    if(nWarnings==maxWarnings) {
      G4cout << "[INCL++] INCL++ has already emitted " << maxWarnings << " warnings and will emit no more." << G4endl;
    }
  }
}

void G4INCLXXInterfaceStore::EmitBigWarning(const G4String &message) const {
  G4cout
    << G4endl
    << "================================================================================"
    << G4endl
    << "                                 INCL++ WARNING                                 "
    << G4endl
    << message
    << G4endl
    << "================================================================================"
    << G4endl
    << G4endl;
}

