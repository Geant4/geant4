// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3PartTable.hh,v 1.1 1999-01-07 16:06:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/gdlist.h>

#include "globals.hh"
#include "G4ParticleDefinition.hh"

struct PartTableEntry {
    G4int partid;
    G4ParticleDefinition* partpt;
};
declare (RWGDlist, PartTableEntry)

class G3PartTable {
private:
    RWGDlist(PartTableEntry) PartT;
    RWGDlist(PartTableEntry)* PartTable;
public:
    G3PartTable();
    ~G3PartTable();
    G4ParticleDefinition* get(G4int partid);
    void put(G4int* partid, G4ParticleDefinition* partpt);
};

extern G3PartTable G3Part;
