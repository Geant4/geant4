// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTable.cc,v 1.1 1999-01-07 16:06:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3DetTable.hh"

class G4VSensitiveDetector;

RWBoolean DetTableMatch(const DetTableEntry *DetTentry, const void *pt)
{
    G4String *name;
    name = (G4String*) pt;
    if (DetTentry->name == *name)
        {
            return TRUE;
        } else {
            return FALSE;
        }
}

G3DetTable::G3DetTable()
{
    DetTable = &DetT;
}

G3DetTable::~G3DetTable()
{
    while (! DetTable->isEmpty()) {
        DetTableEntry *DetTentry = DetTable->last();
        DetTable->removeReference(DetTentry);
        delete DetTentry;
    }
    delete DetTable;
}

G4VSensitiveDetector *G3DetTable::get(G4String detname)
{
    const void *pt;
    pt = &detname;
    DetTableEntry *DetTentry = DetTable->find(DetTableMatch, pt);
    if (DetTentry == NULL) {
        return NULL;
    } else {
        return DetTentry->detpt;
    }
}

G4int G3DetTable::GetID(G4String detname)
{
    const void *pt;
    pt = &detname;
    DetTableEntry *DetTentry = DetTable->find(DetTableMatch, pt);
    if (DetTentry == NULL) {
        return 0;
    } else {
        return DetTentry->detid;
    }
}

void G3DetTable::put(G4String name, G4int detid, G4VSensitiveDetector *detpt)
{
    DetTableEntry *DetTentry = new DetTableEntry;
    DetTentry->name = name;
    DetTentry->detid = detid;
    DetTentry->detpt = detpt;
    DetTable->append(DetTentry);
}
