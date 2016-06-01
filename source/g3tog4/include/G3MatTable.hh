// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.hh,v 2.1 1998/07/12 02:54:12 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/gdlist.h>

#include "globals.hh"
#include "G4Material.hh"

struct MatTableEntry {
    G4int matid;
    G4Material* matpt;
};
declare (RWGDlist, MatTableEntry)

class G3MatTable {
private:
    RWGDlist(MatTableEntry) MatT;
    RWGDlist(MatTableEntry)* MatTable;
public:
    G3MatTable();
    ~G3MatTable();
    G4Material* get(G4int matid);
    void put(G4int* matid, G4Material* matpt);
};

extern G3MatTable G3Mat;
