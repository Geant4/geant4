// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.cc,v 1.10 1999-12-05 17:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MedTable.hh"

G3MedTable::G3MedTable()
{
  fMedVector = new G3MediumVector();
}

G3MedTable::~G3MedTable()
{
  fMedVector->clearAndDestroy();
  delete fMedVector;
}

G3MedTableEntry* G3MedTable::get(G4int id) const
{
  for (G4int i=0; i< fMedVector->entries(); i++) {
    G3MedTableEntry* mte = (*fMedVector)[i];
    if (id == mte->GetID()) return mte;
  }
  return 0;
}    

void G3MedTable::put(G4int id, G4Material* material, G4MagneticField* field,
       G4UserLimits* limits, G4int isvol)
{
  G3MedTableEntry* mte 
    = new G3MedTableEntry(id, material, field, limits, isvol);
  fMedVector->insert(mte);
}

void G3MedTable::Clear()
{
  fMedVector->clearAndDestroy();
}
