// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.1 1999-01-07 16:06:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G3RotTable.hh"

RWBoolean RotTableMatch(const RotTableEntry *RotTentry, const void *pt)
{
    G4int rotid;
    rotid = *((G4int*) pt);
    if (RotTentry->rotid == rotid)
        {
            return TRUE;
        } else {
            return FALSE;
        }
}

G3RotTable::G3RotTable()
{
    RotTable = &RotT;
}

G3RotTable::~G3RotTable()
{
    while (! RotTable->isEmpty()) {
        RotTableEntry *RotTentry = RotTable->last();
        RotTable->removeReference(RotTentry);
        delete RotTentry;
    }
    delete RotTable;
}

G4RotationMatrix *G3RotTable::get(G4int rotid)
{
    const void *pt;
    if ( rotid == 0 ) return NULL;
    pt = &rotid;
    RotTableEntry *RotTentry = RotTable->find(RotTableMatch, pt);
    if (RotTentry == NULL) {
        return NULL;
    } else {
        return RotTentry->rotpt;
    }
}

void G3RotTable::put(G4int *rotid, G4RotationMatrix *rotpt)
{
    RotTableEntry *RotTentry = new RotTableEntry;
    RotTentry->rotid = *rotid;
    RotTentry->rotpt = rotpt;
    RotTable->append(RotTentry);
}
