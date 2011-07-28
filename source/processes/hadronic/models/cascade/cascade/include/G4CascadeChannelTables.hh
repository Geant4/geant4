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
// Factory function to return pointer to Bertini cross-section table based on
// collision initial state (hadron type codes).
//
// Author:  Michael Kelsey (SLAC)

#ifndef G4_CASCADE_CHANNEL_TABLES_HH
#define G4_CASCADE_CHANNEL_TABLES_HH

#include "globals.hh"

class G4CascadeChannel;


namespace G4CascadeChannelTables {
  // Argument is interaction code, product of G4InuclEP types
  const G4CascadeChannel* GetTable(G4int initialState);

  // Arguments are individual G4InuclElementaryParticle types
  inline const G4CascadeChannel* GetTable(G4int had1, G4int had2);

  // Convenience function for diagnostic output
  void PrintTable(G4int initialState);
}

inline const G4CascadeChannel* 
G4CascadeChannelTables::GetTable(G4int had1, G4int had2) {
  return GetTable(had1*had2);
}

#endif	/* G4_CASCADE_CHANNEL_TABLES_HH */
