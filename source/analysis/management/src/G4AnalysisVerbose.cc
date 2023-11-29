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

// Author: Ivana Hrivnacova, 17/10/2011  (ivana@ipno.in2p3.fr)

#include "G4AnalysisVerbose.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

//_____________________________________________________________________________
G4AnalysisVerbose::G4AnalysisVerbose() = default;

//
// public method
//

//_____________________________________________________________________________
void G4AnalysisVerbose::Message(G4int level,
                                const G4String& action,
                                const G4String& object,
                                const G4String& objectName,
                                G4bool success) const
{
  if ( level == 0 ) return;

  if ( level < 0 || level > fkMaxLevel ) {
     // add exception
     return;
  }

  G4cout << "... "
         << fToBeDoneText[level-1]
         << action
         << " "
         << object;
  if (objectName.size() != 0u) {
     G4cout << " : " << objectName;
  }

  if (success) {
     G4cout << " " << fDoneText[level-1];
  }
  else {
     G4cout << " " << fFailureText;
  }

  G4cout << G4endl;
}
