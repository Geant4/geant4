// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.cc,v 1.1 1999-01-07 16:06:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3MedTable.hh"

RWBoolean MedTableMatch(const MedTableEntry *MedTentry, const void *pt)
{
    G4int medid;
    medid = *((G4int*) pt);
    if (MedTentry->medid == medid)
        {
            return TRUE;
        } else {
            return FALSE;
        }
}

G3MedTable::G3MedTable()
{
    MedTable = &MedT;
}

G3MedTable::~G3MedTable()
{
    while (! MedTable->isEmpty()) {
        MedTableEntry *MedTentry = MedTable->last();
        MedTable->removeReference(MedTentry);
        delete MedTentry;
    }
    delete MedTable;
}

G4Material *G3MedTable::GetMat(G4int medid)
{
    const void *pt;
    pt = &medid;
    MedTableEntry *MedTentry = MedTable->find(MedTableMatch, pt);
    if (MedTentry == NULL) {
        return NULL;
    } else {
        return MedTentry->matpt;
    }
}

G4MagneticField *G3MedTable::GetMag(G4int medid)
{
    const void *pt;
    pt = &medid;
    MedTableEntry *MedTentry = MedTable->find(MedTableMatch, pt);
    if (MedTentry == NULL) {
        return NULL;
    } else {
        return MedTentry->magpt;
    }
}

G4UserLimits *G3MedTable::GetLim(G4int medid)
{
    const void *pt;
    pt = &medid;
    MedTableEntry *MedTentry = MedTable->find(MedTableMatch, pt);
    if (MedTentry == NULL) {
        return NULL;
    } else {
        return MedTentry->limpt;
    }
}

void G3MedTable::put( G4int medid, G4Material *matpt, G4MagneticField *magpt,
                      G4UserLimits *limpt, G4int isvol, 
                      G4double deemax, G4double epsil)
{
    MedTableEntry *MedTentry = new MedTableEntry;
    MedTentry->medid = medid;
    MedTentry->matpt = matpt;
    MedTentry->magpt = magpt;
    MedTentry->limpt = limpt;
    MedTentry->isvol = isvol;
    MedTentry->deemax = deemax;
    MedTentry->epsil = epsil;
    MedTable->insert(MedTentry);
}
