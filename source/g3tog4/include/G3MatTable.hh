#ifndef _G3MATTABLE_
#define _G3MATTABLE_ 1

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.hh,v 1.3 1999-05-06 17:46:37 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.

#include <rw/tpordvec.h>

class G4Material;

class G3MatTable {
private:
  RWTPtrOrderedVector<G4Material>* _Mat;
public:
  G3MatTable();
  ~G3MatTable();
  G4Material* get(G4int matid);
  void put(G4int matid, G4Material* matpt);
};

extern G3MatTable G3Mat;
#endif
