// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTable.hh,v 1.1 1999-01-07 16:06:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/gdlist.h>

#include "globals.hh"

class G4VSensitiveDetector;

struct DetTableEntry {
    G4String name;
    G4int detid;
    G4VSensitiveDetector* detpt;
};
declare (RWGDlist, DetTableEntry)

class G3DetTable {
private:
    RWGDlist(DetTableEntry) DetT;
    RWGDlist(DetTableEntry)* DetTable;
public:
    G3DetTable();
    ~G3DetTable();
    G4VSensitiveDetector* get(G4String name);
    G4int GetID(G4String name);
    void put(G4String name, G4int detid, G4VSensitiveDetector* detpt);
};

extern G3DetTable G3Det;
