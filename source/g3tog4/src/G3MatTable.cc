// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.1 1999-01-07 16:06:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3MatTable.hh"

RWBoolean MatTableMatch(const MatTableEntry* MatTentry, const void* pt)
{
    G4int matid;
    matid = *((G4int*) pt);
    if (MatTentry->matid == matid)
        {
            return TRUE;
        } else {
            return FALSE;
        }
}

G3MatTable::G3MatTable()
{
    MatTable = &MatT;
}

G3MatTable::~G3MatTable()
{
    while (! MatTable->isEmpty()) {
        MatTableEntry* MatTentry = MatTable->last();
        MatTable->removeReference(MatTentry);
        delete MatTentry;
    }
    delete MatTable;
}

G4Material* G3MatTable::get(G4int matid)
{
    const void* pt;
    pt = &matid;
    MatTableEntry* MatTentry = MatTable->find(MatTableMatch, pt);
    if (MatTentry == NULL) {
        return NULL;
    } else {
        return MatTentry->matpt;
    }
}

void G3MatTable::put(G4int* matid, G4Material* matpt)
{
    MatTableEntry* MatTentry = new MatTableEntry;
    MatTentry->matid = *matid;
    MatTentry->matpt = matpt;
    MatTable->append(MatTentry);
}
