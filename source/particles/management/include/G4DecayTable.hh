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
// $Id: G4DecayTable.hh 92057 2015-08-14 13:34:55Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      7 July 1996 H.Kurashige
//
// ----------------------------------------
//      implementation for STL          14 Feb. 2000 H.Kurashige
//      
// ------------------------------------------------------------

#ifndef G4DecayTable_h
#define G4DecayTable_h 1

#include "G4ios.hh"
#include <vector>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDecayChannel.hh"

class G4DecayTable 
{
 // Class Description
 //   G4DecayTable is the table of pointer to G4VDecayChannel.
 //   Decay channels inside is sorted by using decay branching ratio
 //

 public:
  typedef std::vector<G4VDecayChannel*> G4VDecayChannelVector;

  //constructors
 public:
    G4DecayTable();
    ~G4DecayTable();

 private:
  // hide copy constructor and assignment operator by declaring "private"
  //  (Implementation does not make sense )
    G4DecayTable(const G4DecayTable &){};
    G4DecayTable & operator=(const G4DecayTable &){return *this;};

 public:
    // equality operators
    G4int operator==(const G4DecayTable &right) const {return (this == &right);};
    G4int operator!=(const G4DecayTable &right) const {return (this != &right);};

 public: // With Description
    void  Insert( G4VDecayChannel* aChannel);
    // Insert a decay channel at proper position 
    // (i.e. sorted by using branching ratio ) 

    G4int entries() const;
    // Returns number of decay channels inside

 public: // With Description
    G4VDecayChannel* SelectADecayChannel(G4double parentMass= -1.);
    // A decay channel is selected at random according to the branching ratio 
  
    G4VDecayChannel* GetDecayChannel(G4int index) const;
    G4VDecayChannel* operator[](G4int index);
    // Get index-th Decay channel
 
    void DumpInfo() const;

 private:
    G4ParticleDefinition       *parent;
    G4VDecayChannelVector       *channels;
};

inline     
 G4int G4DecayTable::entries() const
{
  return channels->size();
}

inline     
 G4VDecayChannel* G4DecayTable::operator[](G4int index)
{
  return (*channels)[index];
}

 
inline     
 G4VDecayChannel* G4DecayTable::GetDecayChannel(G4int index) const
{
  G4VDecayChannel* selectedChannel = 0;
  if ( (index>=0) && (index<G4int(channels->size())) ){
    selectedChannel = (*channels)[index];
  }
  return selectedChannel;
}
 
 
#endif
