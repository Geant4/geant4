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
//
// by I.Hrivnacova, 27 Sep 99

#include "G3MedTable.hh"

G3MedTable::G3MedTable()
{
  fMedVector = new G3MediumVector();
}

G3MedTable::~G3MedTable()
{
  Clear();
  delete fMedVector;
}

G3MedTableEntry* G3MedTable::get(G4int id) const
{
  for (size_t i=0; i< fMedVector->size(); i++) {
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
  fMedVector->push_back(mte);
}

G3MedTableEntry* G3MedTable::GetMTE(G4int i) const
{
  if (i<0 || i>= G4int(fMedVector->size())) 
    return 0;
  
  return (*fMedVector)[i];
}    
    
G4int G3MedTable::GetSize() const
{
  return (G4int)fMedVector->size();
}    

void G3MedTable::Clear()
{
  G3MedTableEntry* a;
  while (fMedVector->size()>0) {
    a = fMedVector->back();
    fMedVector->pop_back();
    for (G3MediumVector::iterator i=fMedVector->begin();
                                  i!=fMedVector->end();){
      if (*i==a) {
	i = fMedVector->erase(i);
      }
      else {
	++i;
      }
    } 
    if ( a )  delete a;    
  } 
}
