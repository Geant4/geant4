// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTableEntry.cc,v 1.2 1999-12-05 17:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MedTableEntry.hh"

#include "G4Material.hh"
#include "G4MagneticField.hh"
#include "G4UserLimits.hh"

G3MedTableEntry::G3MedTableEntry(G4int id, G4Material* material, 
       G4MagneticField* field, G4UserLimits* limits, G4int isvol)
  : fID(id),
    fMaterial(material),
    fField(field),
    fLimits(limits),
    fISVOL(isvol)
{}

G3MedTableEntry::G3MedTableEntry(const G3MedTableEntry& right)
  : fID(right.GetID()),
    fMaterial(right.GetMaterial()),
    fField(right.GetField()),
    fLimits(right.GetLimits()),
    fISVOL(right.GetISVOL())    
{}    

G3MedTableEntry::~G3MedTableEntry()
{}

const G3MedTableEntry& 
G3MedTableEntry::operator=(const G3MedTableEntry& right)
{ 
  fID = right.GetID();
  fMaterial = right.GetMaterial();     
  fField = right.GetField();
  fLimits = right.GetLimits();
  fISVOL = right.GetISVOL();  
  return *this;
}

G4int G3MedTableEntry::operator==(const G3MedTableEntry& right) const
{ 
  if (fID == right.GetID()) 
    return 1;
  else
    return 0;
}

G4int G3MedTableEntry::operator!=(const G3MedTableEntry& right) const
{ 
  if (*this == right) 
    return 0;
  else
    return 1;
}

