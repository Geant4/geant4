// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.hh,v 1.6 1999-12-05 17:50:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3_MED_TABLE_H
#define G3_MED_TABLE_H

#include "G3MedTableEntry.hh"

#include "globals.hh"

#include "g4rw/tpordvec.h"

class G4Material;
class G4MagneticField;
class G4UserLimits;

typedef G4RWTPtrOrderedVector<G3MedTableEntry>  G3MediumVector;

class G3MedTable
{
  public:
    G3MedTable();
    virtual ~G3MedTable();
    
    // methods
    G3MedTableEntry* get(G4int id) const;
    void put(G4int id, G4Material* material, G4MagneticField* field,
       G4UserLimits* limits, G4int isvol);
    void Clear();

  private:
    G3MediumVector*  fMedVector;
};

extern G3MedTable G3Med;

#endif //G3_MED_TABLE_H
