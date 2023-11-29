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

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#include "G4AnalysisManagerState.hh"

#include "G4AnalysisUtilities.hh"
#include "G4UnitsTable.hh"

#include <utility>

using namespace G4Analysis;

//_____________________________________________________________________________
G4AnalysisManagerState::G4AnalysisManagerState(G4String type, G4bool isMaster)
  : fType(std::move(type)),
    fIsMaster(isMaster),
    fThreadId(G4Threading::G4GetThreadId())
{}

//
// private methods
//

//_____________________________________________________________________________
void G4AnalysisManagerState::SetVerboseLevel(G4int verboseLevel)
{
  if ( verboseLevel == fVerboseLevel ) return;

  if ( verboseLevel < 0 ) {
    Warn("Cannot set value < 0", fkClass, "SetVerboseLevel");
    return;
  }

  fVerboseLevel = verboseLevel;
}

//
// public methods
//

//_____________________________________________________________________________
void G4AnalysisManagerState::Message(
  [[maybe_unused]] G4int level,
  [[maybe_unused]] const G4String& action,
  [[maybe_unused]] const G4String& objectType,
  [[maybe_unused]] const G4String& objectName,
  [[maybe_unused]] G4bool success ) const
{
#ifdef G4VERBOSE
  // Skip message if of higher level than that is set
  if (fVerboseLevel < level) return;

  // Print message
  fVerbose.Message(level, action, objectType, objectName, success);
#endif
}
