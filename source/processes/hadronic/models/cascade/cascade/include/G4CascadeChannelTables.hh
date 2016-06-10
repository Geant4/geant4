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
// $Id: G4CascadeChannelTables.hh 69336 2013-04-30 20:20:23Z mkelsey $
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
// 20130129  M. Kelsey -- Drop load-on-demand interfaces, fill in ctor
// 20130306  M. Kelsey -- Add inclusive printing of all tables here
// 20130429  M. Kelsey -- Change instance to thread-local pointer.
// 20141121  M. Kelsey -- Dtor must be public for end-of-job cleanup.

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

  // Convenience functions for diagnostic output
  static void Print(std::ostream& os=G4cout);
  static void PrintTable(G4int initialState, std::ostream& os=G4cout);

public:
  ~G4CascadeChannelTables();

private:
  static const G4CascadeChannelTables& instance();		// Singleton
  static G4ThreadLocal G4CascadeChannelTables* theInstance;	// per thread

  G4CascadeChannelTables();

  // Fetch table from map if already registered, or return null
  const G4CascadeChannel* FindTable(G4int initialState) const;

  typedef std::map<G4int, G4CascadeChannel*> TableMap;
  TableMap tables;
};

#endif	/* G4_CASCADE_CHANNEL_TABLES_HH */
