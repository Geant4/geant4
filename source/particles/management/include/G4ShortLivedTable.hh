//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ShortLivedTable.hh,v 1.10 2001-07-11 10:01:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added Remove()                  06 Nov.,98 H.Kurashige

#ifndef G4ShortLivedTable_h
#define G4ShortLivedTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "g4std/vector"

class G4ParticleTable;

class G4ShortLivedTable
{
 // Class Description
 //   G4ShortLivedTable is the table of pointer to G4ParticleDefinition
 //   In G4ShortLivedTable, each G4ParticleDefinition pointer is stored
 //

 public:
   // Use STL Vector as list of shortlives
   typedef G4std::vector<G4ParticleDefinition*>  G4ShortLivedList;

 public:
   G4ShortLivedTable();

 protected:
   G4ShortLivedTable(const  G4ShortLivedTable &right);

 public: // With Description
   virtual ~G4ShortLivedTable();

   G4bool                IsShortLived(G4ParticleDefinition*) const;
   // return true if the particle is shortlived particle
  
   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 

   G4int                 Entries() const;
   G4bool                Contains(const G4ParticleDefinition *particle) const;
   void                  Insert(G4ParticleDefinition* particle);
   void                  Remove(G4ParticleDefinition* particle);
   G4ParticleDefinition* GetParticle(G4int index) const;

 protected:
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


inline 
 G4ParticleDefinition*  G4ShortLivedTable::GetParticle(G4int index) const
{
  if ( (index >=0 ) && (index < Entries()) ) {
    return (*fShortLivedList)[index];
  } else {
    return 0; 
  } 
}



#endif
