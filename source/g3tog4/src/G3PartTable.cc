// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3PartTable.cc,v 1.1 1999-01-07 16:06:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3PartTable.hh"

RWBoolean PartTableMatch(const PartTableEntry *PartTentry, const void *pt)
{
    G4int partid;
    partid = *((G4int*) pt);
    if (PartTentry->partid == partid)
        {
            return TRUE;
        } else {
            return FALSE;
        }
}

G3PartTable::G3PartTable()
{
    PartTable = &PartT;
}

G3PartTable::~G3PartTable()
{
    while (! PartTable->isEmpty()) {
        PartTableEntry *PartTentry = PartTable->last();
        PartTable->removeReference(PartTentry);
        delete PartTentry;
    }
    delete PartTable;
}

G4ParticleDefinition *G3PartTable::get(G4int partid)
{
    const void *pt;
    pt = &partid;
    PartTableEntry *PartTentry = PartTable->find(PartTableMatch, pt);
    if (PartTentry == NULL) {
        return NULL;
    } else {
        return PartTentry->partpt;
    }
}

void G3PartTable::put(G4int *partid, G4ParticleDefinition *partpt)
{
    PartTableEntry *PartTentry = new PartTableEntry;
    PartTentry->partid = *partid;
    PartTentry->partpt = partpt;
    PartTable->append(PartTentry);
}
