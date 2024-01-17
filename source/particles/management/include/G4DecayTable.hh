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
// G4DecayTable
//
// Class description:
//
// G4DecayTable is the table of pointers to G4VDecayChannel.
// Decay channels inside are sorted by using decay branching ratio

// Author: H.Kurashige, 7 July 1996
// --------------------------------------------------------------------
#ifndef G4DecayTable_hh
#define G4DecayTable_hh 1

#include "G4ParticleDefinition.hh"
#include "G4VDecayChannel.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <vector>

class G4DecayTable
{
  public:
    using G4VDecayChannelVector = std::vector<G4VDecayChannel*>;

    G4DecayTable();
    ~G4DecayTable();

    G4DecayTable(const G4DecayTable&) = delete;
    G4DecayTable& operator=(const G4DecayTable&) = delete;

    // Equality operators
    inline G4bool operator==(const G4DecayTable& right) const;
    inline G4bool operator!=(const G4DecayTable& right) const;

    // Insert a decay channel at proper position
    // (i.e. sorted by using branching ratio )
    void Insert(G4VDecayChannel* aChannel);

    // Returns number of decay channels inside
    inline G4int entries() const;

    // A decay channel is selected at random according to the branching ratio
    G4VDecayChannel* SelectADecayChannel(G4double parentMass = -1.);

    // Get index-th decay channel
    inline G4VDecayChannel* GetDecayChannel(G4int index) const;
    inline G4VDecayChannel* operator[](G4int index);

    void DumpInfo() const;

  private:
    G4ParticleDefinition* parent = nullptr;
    G4VDecayChannelVector* channels = nullptr;
};

// ------------------------
// Inline methods
// ------------------------

inline G4bool G4DecayTable::operator==(const G4DecayTable& right) const
{
  return (this == &right);
}

inline G4bool G4DecayTable::operator!=(const G4DecayTable& right) const
{
  return (this != &right);
}

inline G4int G4DecayTable::entries() const
{
  return G4int(channels->size());
}

inline G4VDecayChannel* G4DecayTable::operator[](G4int index)
{
  return (*channels)[index];
}

inline G4VDecayChannel* G4DecayTable::GetDecayChannel(G4int index) const
{
  G4VDecayChannel* selectedChannel = nullptr;
  if ((index >= 0) && (index < G4int(channels->size()))) {
    selectedChannel = (*channels)[index];
  }
  return selectedChannel;
}

#endif
