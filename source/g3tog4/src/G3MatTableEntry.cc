// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTableEntry.cc,v 1.2 1999-12-05 17:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MatTableEntry.hh"

#include "G4Material.hh"

G3MatTableEntry::G3MatTableEntry(G4int id, G4Material* material)
  : fID(id),
    fMaterial(material)
{}

G3MatTableEntry::G3MatTableEntry(const G3MatTableEntry& right)
  : fID(right.GetID()),
    fMaterial(right.GetMaterial())
{}    

G3MatTableEntry::~G3MatTableEntry()
{}

const G3MatTableEntry& 
G3MatTableEntry::operator=(const G3MatTableEntry& right)
{ 
  fID = right.GetID();
  fMaterial = right.GetMaterial();     
  return *this;
}

G4int G3MatTableEntry::operator==(const G3MatTableEntry& right) const
{ 
  if (fID == right.GetID()) 
    return 1;
  else
    return 0;
}

G4int G3MatTableEntry::operator!=(const G3MatTableEntry& right) const
{ 
  if (*this == right) 
    return 0;
  else
    return 1;
}

