// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.hh,v 1.1 1999-01-07 16:06:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 logical volume table.
// Maps G3 logical volume names to their G4 logical volume object
//   counterparts.
// Maintains a linked List of G3 volume names//G4 logical volume pointer
//   pairs.

#include <rw/gdlist.h>
#include <rw/tpordvec.h>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"

struct VolTableEntry {
    G4String vname;
    G4VSolid* solid;
    G4LogicalVolume* lvpt;
    G4VPhysicalVolume* pvpt;
    G4double rangehi[3];   // ranges 
    G4double rangelo[3];   // ranges 
    EAxis *axiscode;       // axis codes
    G4int tmed;
    G4String shape;
    G4double* par;
    G4int npar;
    G4int count; // number of instantiations for vols indexed due to -ve pars
    RWTPtrOrderedVector<G4LogicalVolume> lVols; // vector of instantiated volumes
    RWTPtrOrderedVector<G4VPhysicalVolume> pVols; // vector of instantiated phys volumes
};
declare (RWGDlist, VolTableEntry)

class G3VolTable {
private:
    RWGDlist(VolTableEntry) VolT;
    RWGDlist(VolTableEntry)* VolTable;
    VolTableEntry* curEntry;
    G4int ScanTmed;
    G4int nEntry;
    EAxis Rect[3];
    EAxis Polar[3];
    static G4VPhysicalVolume* mothPV; // mother of all mothers
    static G4LogicalVolume* mothLV; // mother of all mothers logical volume
    RWTPtrOrderedVector<G4VPhysicalVolume>* pVolsPtr;
public:
    G3VolTable();
    ~G3VolTable();
    G4LogicalVolume* GetLV();
    G4LogicalVolume* GetLVx(G4String &);
    G4LogicalVolume* GetLV(G4String* vname);
    G4VPhysicalVolume* GetPV(); // simply return the mother of all volumes
    G4VPhysicalVolume* GetPV(G4int npv); // Use pVols table of last retrieved phys vol
        // FindPV : look up PV in table to use in placement; if no associated PV, use 
        //          global mother
    void FindPV(G4int* npv, G4String* vname);
        // MatchPV: look up PV only
    void MatchPV(G4int* npv, G4String* vname);
    void GetLVPars(G4String* vname, G4int iaxis, G4double* rangehi, 
                   G4double* rangelo, EAxis*,
                   G4String* shape, G4int* nmed, G4double* par[],
                   G4int* npar, G4VSolid* solid);
    void GetLVInfo(G4String* vname, G4String* shape, G4int* nmed);
    void AddConstituentLVol(G4String* vname, G4LogicalVolume* lvol);
    void PutLV(G4String* vname, G4LogicalVolume* lvpt, G4int tmed,
               G4double rnghi[], G4double rnglo[],
               EAxis axis, G4String shape,
               G4double par[], G4int npar, G4VSolid* solid);
         // For case where nothing is known about the volume; it doesn't
         // come from g3tog4 (eg. externally specified global mother)
    void PutLV(G4String* vname, G4LogicalVolume* lvpt);
        // Builds List of all physical volume instantiations of a Indexed volume implemented
        // as individual logical volumes (so when daughters are added, they must be added
        // to the logical volumes of all the physical volumes in this List)
    void PutPV(G4String* vname, G4VPhysicalVolume* lvpt);
        // For non-Indexed volumes, builds a List of one element only
    void PutPV1(G4String* vname, G4VPhysicalVolume* lvpt);
    void MatchTmed(G4int tmed);
    G4LogicalVolume* NextTmed();           
    G4int GetEntryCount() { return VolTable->entries();};
    void SetMother(G4VPhysicalVolume* mpv); // set mother-of-all
    void SetMother(G4LogicalVolume* mlv); // set mother-of-all
    void SetMother(); // set mother of all pointers to null
    VolTableEntry* find(G4String* vname); // find entry from vol name
    G4int Increment(G4String* vname); // Increment instantiation counter
};

extern G3VolTable G3Vol; // logical volumes


