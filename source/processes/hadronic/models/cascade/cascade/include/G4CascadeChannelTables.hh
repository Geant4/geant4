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
// Factory class to return pointer to Bertini cross-section table based on
// collision initial state (hadron type codes).
//
// Author:  Michael Kelsey (SLAC)
//
// 20110728  M. Kelsey -- Use static instance() function to work around
//		"disappearance" bug on Linux (GCC 4.1.2).
// 20110915  M. Kelsey -- Move static implementations (won't work on Windows);
//		Add local table instantiating function (to replace self-reg)
// 20110923  M. Kelsey -- Add optional stream& argument to printTable().

#ifndef G4_CASCADE_CHANNEL_TABLES_HH
#define G4_CASCADE_CHANNEL_TABLES_HH

#include "globals.hh"
#include <iosfwd>
#include <map>

class G4CascadeChannel;


class G4CascadeChannelTables {
public:
  // Argument is interaction code, product of G4InuclEP types
  static const G4CascadeChannel* GetTable(G4int initialState);

  // Arguments are individual G4InuclElementaryParticle types
  static const G4CascadeChannel* GetTable(G4int had1, G4int had2);

  // Convenience function for diagnostic output
  static void PrintTable(G4int initialState, std::ostream& os=G4cout);

  // Register cross-section table for later lookup
  static void AddTable(G4int initialState, G4CascadeChannel* table);

private:
  static G4CascadeChannelTables& instance();	// Singleton

  G4CascadeChannelTables() { tables.clear(); }	// Must initialize map
  ~G4CascadeChannelTables() {}	// Tables are created externally, not owned

  // Fetch table from map if already registered, or return null
  const G4CascadeChannel* FindTable(G4int initialState);

  // Save table for specified interaction in map 
  void SaveTable(G4int initialState, G4CascadeChannel* table);

  // Special function to create and store table for specified interaction
  const G4CascadeChannel* LoadTable(G4int initialState);

  typedef std::map<G4int, G4CascadeChannel*> TableMap;
  typedef std::pair<const G4int, G4CascadeChannel*> TableEntry;
  TableMap tables;

  static void DeleteTable(TableEntry& t);	// For use by destructor
};

#endif	/* G4_CASCADE_CHANNEL_TABLES_HH */
