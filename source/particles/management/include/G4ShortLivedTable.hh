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
// $Id: G4ShortLivedTable.hh,v 1.14 2010-08-10 15:47:42 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added clear()                   20 Mar.,08 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige

#ifndef G4ShortLivedTable_h
#define G4ShortLivedTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include <vector>

class G4ParticleTable;

class G4ShortLivedTable
{
 // Class Description
 //   G4ShortLivedTable is the table of pointer to G4ParticleDefinition
 //   In G4ShortLivedTable, each G4ParticleDefinition pointer is stored
 //

 public:
   // Use STL Vector as list of shortlives
   typedef std::vector<const G4ParticleDefinition*>  G4ShortLivedList;

 public:
   G4ShortLivedTable();

 protected:
   G4ShortLivedTable(const  G4ShortLivedTable &right);

 public: // With Description
   virtual ~G4ShortLivedTable();

   G4bool                IsShortLived(const G4ParticleDefinition*) const;
   // return true if the particle is shortlived particle
  
   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 

   G4int                 Entries() const;
   // return number of particles in the list

   G4bool                Contains(const G4ParticleDefinition *particle) const;
   // return true if the list contains the specified particle 

   void                  Insert(const G4ParticleDefinition* particle);
   // add the particle in the list

   void                  Remove(const G4ParticleDefinition* particle);
   // remove the particle (not delete) from the list 

   G4ParticleDefinition* GetParticle(G4int index) const;
   // return the i-th particle in the list

   G4int                 size() const;
   // return number of particles in the list
  
   void                  clear();
   // remove all particles (not delete) from the list

 protected://Without Description
   G4int                GetVerboseLevel() const;

 private:
   G4ShortLivedList*                  fShortLivedList;

};

inline G4bool  G4ShortLivedTable::Contains(const G4ParticleDefinition* particle) const
{
  G4ShortLivedList::iterator i;
  for (i = fShortLivedList->begin(); i!= fShortLivedList->end(); ++i) {
    if (**i==*particle) return true;
  }
  return false;
}

inline G4int G4ShortLivedTable::Entries() const
{
  return fShortLivedList->size();
}

inline G4int G4ShortLivedTable::size() const
{
  return fShortLivedList->size();
}

inline void G4ShortLivedTable::clear()
{
  fShortLivedList->clear();
}

inline 
 G4ParticleDefinition*  G4ShortLivedTable::GetParticle(G4int index) const
{
  if ( (index >=0 ) && (index < Entries()) ) {
    return const_cast<G4ParticleDefinition*>( (*fShortLivedList)[index] );
  } else {
    return 0; 
  } 
}



#endif
