// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.12 1999-12-05 17:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MatTable.hh"

G3MatTable::G3MatTable()
{
  fMatVector = new G3MaterialVector();
}

G3MatTable::~G3MatTable()
{
  fMatVector->clearAndDestroy();
  delete fMatVector;
}

G4Material* G3MatTable::get(G4int id) const
{
  for (G4int i=0; i< fMatVector->entries(); i++) {
    G3MatTableEntry* mte = (*fMatVector)[i];
    if (id == mte->GetID()) return mte->GetMaterial();
  }
  return 0;
}    

void G3MatTable::put(G4int id, G4Material* material)
{
  G3MatTableEntry* mte = new G3MatTableEntry(id, material);
  fMatVector->insert(mte);
}

void G3MatTable::Clear()
{
  fMatVector->clearAndDestroy();
}
