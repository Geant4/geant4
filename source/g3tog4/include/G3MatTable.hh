// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.hh,v 1.9 2000-11-24 09:50:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 materials table.
// Maps G3 material indices to their G4Material object counterparts.

// ----------------------
//
// by I.Hrivnacova, 27 Sep 99

#ifndef G3MATTABLE_HH
#define G3MATTABLE_HH 1

#include "G3MatTableEntry.hh"

#include "globals.hh"

#include "g4rw/tpordvec.h"

class G4Material;

typedef G4RWTPtrOrderedVector<G3MatTableEntry>  G3MaterialVector;

class G3MatTable
{
  public: // with description

    G3MatTable();
    virtual ~G3MatTable();
    
    // methods
    G4Material* get(G4int id) const;
    void put(G4int id, G4Material* material);
    void Clear();

  private:

    G3MaterialVector*  fMatVector;
};

extern G3MatTable G3Mat;

#endif
