// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.13 2001-05-30 08:13:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3RotTable.hh"

G3RotTable::G3RotTable()
{
  fRotVector = new G3RotMatrixVector();
}

G3RotTable::~G3RotTable()
{
  fRotVector->clear();
  //fRotVector->clearAndDestroy();
  delete fRotVector;
}

G4RotationMatrix* G3RotTable::Get(G4int id) const
{
  for (size_t i=0; i<fRotVector->entries(); i++) {
    G3RotTableEntry* rte = (*fRotVector)[i];
    if (id == rte->GetID()) return rte->GetMatrix();
  }
  return 0;
}    

void G3RotTable::Put(G4int id, G4RotationMatrix* matrix)
{
  G3RotTableEntry* rte = new G3RotTableEntry(id, matrix);
  fRotVector->insert(rte);
}

void G3RotTable::Clear()
{
  fRotVector->clearAndDestroy();
}
