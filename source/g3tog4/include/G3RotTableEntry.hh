// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTableEntry.hh,v 1.3 1999-12-09 01:27:46 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3ROTTABLEENTRY_HH
#define G3ROTTABLEENTRY_HH 1

#include "globals.hh"
#include "G4RotationMatrix.hh"

class G3toG4RotationMatrix;

class G3RotTableEntry 
{
  public:
    G3RotTableEntry(G4int id, G4RotationMatrix* matrix);
    G3RotTableEntry(const G3RotTableEntry& right);
    virtual ~G3RotTableEntry();
    
    // operators
    const G3RotTableEntry& operator=(const G3RotTableEntry& right);
    G4int operator==(const G3RotTableEntry& right) const;
    G4int operator!=(const G3RotTableEntry& right) const;

    // get methods
    G4int       GetID() const;
    G4RotationMatrix* GetMatrix() const;
    
  private:
    // data members  
    G4int              fID;
    G4RotationMatrix*  fMatrix;
};

// inline methods

inline G4int G3RotTableEntry::GetID() const
{ return fID; }

inline G4RotationMatrix* G3RotTableEntry::GetMatrix() const
{ return fMatrix; }

#endif
