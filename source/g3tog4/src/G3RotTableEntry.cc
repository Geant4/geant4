// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTableEntry.cc,v 1.2 1999-12-05 17:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3RotTableEntry.hh"

#include "G3toG4RotationMatrix.hh"

G3RotTableEntry::G3RotTableEntry(G4int id, G4RotationMatrix* matrix)
  : fID(id),
    fMatrix(matrix)
{}

G3RotTableEntry::G3RotTableEntry(const G3RotTableEntry& right)
  : fID(right.GetID()),
    fMatrix(right.GetMatrix())
{}    

G3RotTableEntry::~G3RotTableEntry()
{}

const G3RotTableEntry& 
G3RotTableEntry::operator=(const G3RotTableEntry& right)
{ 
  fID = right.GetID();
  fMatrix = right.GetMatrix();     
  return *this;
}

G4int G3RotTableEntry::operator==(const G3RotTableEntry& right) const
{ 
  if (fID == right.GetID()) 
    return 1;
  else
    return 0;
}

G4int G3RotTableEntry::operator!=(const G3RotTableEntry& right) const
{ 
  if (*this == right) 
    return 0;
  else
    return 1;
}

