//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G3MatTable.cc 67982 2013-03-13 10:36:03Z gcosmo $
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MatTable.hh"

G3MatTable::G3MatTable()
{
  fMatVector = new G3MaterialVector();
}

G3MatTable::~G3MatTable()
{
  Clear();
  delete fMatVector;
}

G4Material* G3MatTable::get(G4int id) const
{
  for (size_t i=0; i< fMatVector->size(); i++) {
    G3MatTableEntry* mte = (*fMatVector)[i];
    if (id == mte->GetID()) return mte->GetMaterial();
  }
  return 0;
}    

void G3MatTable::put(G4int id, G4Material* material)
{
  G3MatTableEntry* mte = new G3MatTableEntry(id, material);
  fMatVector->push_back(mte);
}

void G3MatTable::Clear()
{
  G3MatTableEntry* a;
  while (fMatVector->size()>0) {
    a = fMatVector->back();
    fMatVector->pop_back();
    for (G3MaterialVector::iterator i=fMatVector->begin();
                                    i!=fMatVector->end();){
      if (*i==a) {
	i = fMatVector->erase(i);
      }
      else {
	++i;
      }
    } 
    if ( a )  delete a;    
  } 
}
