// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.hh,v 1.1 1999-01-07 16:06:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/gdlist.h>

#include "globals.hh"
#include "G4Material.hh"

class G4MagneticField;
class G4UserLimits;

struct MedTableEntry {
    G4int medid;
    G4Material* matpt;
    G4MagneticField* magpt;
    G4UserLimits* limpt;
    G4int isvol;
    G4double deemax;
    G4double epsil;
};
declare (RWGDlist, MedTableEntry)

class G3MedTable {
private:
    RWGDlist(MedTableEntry) MedT;
    RWGDlist(MedTableEntry)* MedTable;
public:
    G3MedTable();
    ~G3MedTable();
    G4Material* GetMat(G4int medid);  // get associated material
    G4MagneticField* GetMag(G4int medid);
    G4UserLimits* GetLim(G4int medid);
    void put(G4int medid, G4Material* matpt, G4MagneticField* field,
             G4UserLimits* limits, G4int isvol,
             G4double deemax, G4double epsil);
};

extern G3MedTable G3Med;
