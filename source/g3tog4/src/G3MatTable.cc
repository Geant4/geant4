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
// $Id: G3MatTable.cc,v 1.14 2001-07-11 09:58:59 gunter Exp $
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
  for (size_t i=0; i< fMatVector->entries(); i++) {
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
