// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTableEntry.hh,v 1.3 1999-12-09 01:27:43 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MATTABLEENTRY_HH
#define G3MATTABLEENTRY_HH 1

#include "globals.hh"

class G4Material;

class G3MatTableEntry 
{
  public:
    G3MatTableEntry(G4int id, G4Material* material);
    G3MatTableEntry(const G3MatTableEntry& right);
    virtual ~G3MatTableEntry();
    
    // operators
    const G3MatTableEntry& operator=(const G3MatTableEntry& right);
    G4int operator==(const G3MatTableEntry& right) const;
    G4int operator!=(const G3MatTableEntry& right) const;

    // get methods
    G4int       GetID() const;
    G4Material* GetMaterial() const;
    
  private:
    // data members  
    G4int        fID;
    G4Material*  fMaterial;
};

// inline methods

inline G4int G3MatTableEntry::GetID() const
{ return fID; }

inline G4Material* G3MatTableEntry::GetMaterial() const
{ return fMaterial; }

#endif
