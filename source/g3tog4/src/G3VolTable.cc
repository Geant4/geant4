// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTable.cc,v 1.1 1999-01-07 16:06:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3VolTable.hh"
#include "G4LogicalVolumeStore.hh"

G4VPhysicalVolume* G3VolTable::mothPV = NULL; // physical volume of mother of all mothers initialized to NULL.
G4LogicalVolume* G3VolTable::mothLV = NULL; // logical volume of mother of all mothers initialized to NULL.

RWBoolean VolTableMatch(const VolTableEntry *VolTentry, const void *pt)
{
    G4String *vname;
    vname = (G4String*) pt;
    if (VolTentry->vname == *vname)
    {
        return TRUE;
    } else {
        return FALSE;
    }
}

G3VolTable::G3VolTable()
{
    VolTable = &VolT;
    G4int nent = VolTable->entries();
    ScanTmed = 0;
    Rect[0] = kXAxis;
    Rect[1] = kYAxis;
    Rect[2] = kZAxis;
    Polar[0] = kRho;
    Polar[1] = kPhi;
    Polar[2] = kZAxis;
}

G3VolTable::~G3VolTable()
{
    while (! VolTable->isEmpty()) {
        VolTableEntry *VolTentry = VolTable->last();
        VolTable->removeReference(VolTentry);
            // delete VolTentry->par;
        delete VolTentry;
    }
    delete VolTable;
}

G4LogicalVolume* G3VolTable::GetLVx(G4String& vname)
{
        // loop through the G4 List of currently defined logical volumes
    
    RWTPtrOrderedVector <G4LogicalVolume>* _ptrToList;
    _ptrToList = G4LogicalVolumeStore::GetInstance();
    G4LogicalVolume* _lvol;
    G4int nentries = _ptrToList->entries();
    for (G4int i=0; i<nentries; i++){
        _lvol = (*_ptrToList)[i];
        if (_lvol->GetName() == vname) return _lvol;
    }
    return NULL;
}

        
G4LogicalVolume *G3VolTable::GetLV(G4String *vname)
{
    curEntry = find(vname);
    if (curEntry == NULL) {
        return NULL;
    } else {
        return curEntry->lvpt;
    }
}

G4LogicalVolume* G3VolTable::GetLV()
{
//    G3Mother* TV;
//    return TV->GetLV(); // return the pointer to the top-level logical volume
//    G4cout << "LV pointer to top level logical volume" << mothLV << endl;
    return mothLV;
}


G4VPhysicalVolume *G3VolTable::GetPV(G4int npv)
{
    G4int ln = pVolsPtr->length();
    if (npv > ln) {
        G4cout << "G3VolTable::GetPV error: request index too high " << npv <<
            " > " << pVolsPtr->length() << endl;
        return NULL;
    } else {
        if (ln == 0) {
                // No associated physical volume. Use global mother.
            G4cout << "GetPV: returning global mother " << endl;
            return mothPV;
        } else {
            return pVolsPtr->at(npv);
        }
    }
}

G4VPhysicalVolume *G3VolTable::GetPV()
{
    return mothPV;
}

void G3VolTable::FindPV(G4int* nvols, G4String *vname)
{
    const void *pt;
    pt = vname;
    curEntry = find(vname);
    if (curEntry == NULL) {
        *nvols = 0;
        return;
    } else {
        pVolsPtr = &(curEntry->pVols);
        if (pVolsPtr->length() == 0) {
                // use global mother
            pVolsPtr->insert(mothPV);
            G4cout << "Global mother set as physical volume for " << *vname << endl;
        }
        *nvols = pVolsPtr->length();
    }
}

void G3VolTable::MatchPV(G4int* nvols, G4String *vname)
{
    const void *pt;
    pt = vname;
    curEntry = find(vname);
    if (curEntry == NULL) {
        *nvols = -1;
        G4cout << "G3VolTable error: no entry for " << *vname << endl;
        return;
    } else {
        pVolsPtr = &(curEntry->pVols);
        *nvols = pVolsPtr->length();
    }
}

void G3VolTable::GetLVPars(G4String *vname, G4int iaxis, G4double *rnghi,
                           G4double *rnglo,
                           EAxis *axis, G4String *shape, G4int *nmed,
                           G4double *par[], G4int *npar, G4VSolid *solid)
{
    const void *pt;
    pt = vname;
    curEntry = find(vname);
    if (curEntry == NULL) {
        *rnghi = 0.;
        *rnglo = 0.;
            // $$$ *axis = kUndefined;
        *shape = "";
        *nmed = 0;
        *npar = 0;
        return;
    } else {
        *rnghi = curEntry->rangehi[iaxis-1];
        *rnglo = curEntry->rangelo[iaxis-1];
        *axis = curEntry->axiscode[iaxis-1];
        *shape = curEntry->shape;
        *nmed = curEntry->tmed;
        *par = curEntry->par;
        *npar = curEntry->npar;
        solid = curEntry->solid;
    }
}

void G3VolTable::GetLVInfo(G4String *vname, G4String *shape, G4int *nmed)
{
    const void *pt;
    pt = vname;
    curEntry = find(vname);
    if (curEntry == NULL) {
        *shape = "";
        *nmed = 0;
        return;
    } else {
        *shape = curEntry->shape;
        *nmed = curEntry->tmed;
    }
}

void G3VolTable::AddConstituentLVol(G4String* vname, G4LogicalVolume* lvol)
{
    const void *pt;
    pt = vname;
    curEntry = find(vname);
    if (curEntry == NULL) {
        G4cout << "G3 LVol " << *vname << " not found" << endl;
        return;
    } else {
        curEntry->lVols.insert(lvol);
    }
}

void G3VolTable::PutLV(G4String *vname, G4LogicalVolume *lvpt, G4int tmed,
                       G4double rnghi[], G4double rnglo[],
                       EAxis axis, G4String shape,
                       G4double par[], G4int npar, G4VSolid *solid)
{
    G4LogicalVolume *lvptold = NULL;
    if (VolTable->entries() > 0) lvptold = GetLV(vname);
    if (lvptold == NULL) {  // create a new entry
        VolTableEntry *VolTentry = new VolTableEntry;
        VolTentry->vname = *vname;
        VolTentry->lvpt = lvpt;
        VolTentry->tmed = tmed;
        VolTentry->shape = shape;
        VolTentry->npar = npar;
        VolTentry->solid = solid;
        VolTentry->count = 0;
            // associated physical volume is initialized to the root mother
            // physical volume $$$
        VolTentry->pvpt = mothPV;

        G4double *params = new G4double[npar];
        G4int i;
        for (i=0; i<npar; i++) {
            params[i] = par[i];
        }
        VolTentry->par = params;
        if (axis == kXAxis) VolTentry->axiscode = Rect;
        if (axis == kRho) VolTentry->axiscode = Polar;
        for (i=0; i<3; i++) {
            VolTentry->rangehi[i] = rnghi[i];
            VolTentry->rangelo[i] = rnglo[i];
        }
        VolTable->append(VolTentry);
    } else {  // load existing entry
        curEntry->lvpt = lvpt;
    }
}

// For case where nothing is known about the volume; it doesn't
// come from G3toG4 (eg. externally specified global mother)
void G3VolTable::PutLV(G4String *vname, G4LogicalVolume *lvpt)
{
    G4LogicalVolume *lvptold = NULL;
    if (VolTable->entries() > 0) lvptold = GetLV(vname);
    if (lvptold == NULL) {  // create a new entry
        VolTableEntry *VolTentry = new VolTableEntry;
        VolTentry->vname = *vname;
        VolTentry->lvpt = lvpt;
        VolTentry->tmed = 0;
        VolTentry->shape = "unknown";
        VolTentry->npar = 0;
        VolTentry->solid = 0;
        VolTentry->count = 0;
            // associated physical volume is initialized to the root mother
            // physical volume $$$
        VolTentry->pvpt = mothPV;

        VolTentry->par = 0;
        VolTentry->axiscode = Rect;
        G4int i;
        for (i=0; i<3; i++) {
            VolTentry->rangehi[i] = 0.;
            VolTentry->rangelo[i] = 0.;
        }
        VolTable->append(VolTentry);
    } else {  // load existing entry
        curEntry->lvpt = lvpt;
    }
}

void G3VolTable::PutPV(G4String *vname, G4VPhysicalVolume *pvpt)
{
    G4int npv=0;
    MatchPV(&npv,vname);
    if (npv >= 0) {
        curEntry->pVols.insert(pvpt);
    }
}


void G3VolTable::PutPV1(G4String *vname, G4VPhysicalVolume *pvpt)
{
        // Only add it if List is empty
    G4int npv=0;
    MatchPV(&npv,vname);
    if (npv == 0) {
        curEntry->pVols.insert(pvpt);
    }
}

void G3VolTable::MatchTmed(G4int tmed)
{
        // Set tmed to match with NextTmed
    ScanTmed = tmed;
    nEntry = 0;
}

G4LogicalVolume *G3VolTable::NextTmed()
{
        // Return successively all matching lvols; NULL terminates
    for (G4int i=nEntry; i<VolTable->entries(); i++) {
        VolTableEntry *vte = VolTable->at(i);
        if (vte->tmed == ScanTmed) {
            nEntry = nEntry + i + 1;
            return vte->lvpt;
        }
    }
    return NULL;
}

void G3VolTable::SetMother(G4VPhysicalVolume *MotherOfAll)
{
    if (!mothPV) mothPV = MotherOfAll;
}

void G3VolTable::SetMother(G4LogicalVolume *MotherOfAll)
{
    if (!mothLV) mothLV = MotherOfAll;
}

void G3VolTable::SetMother()
{
    mothLV = NULL;
    mothPV = NULL;
}

VolTableEntry *G3VolTable::find(G4String *vname)
{
    G4int nEntries = VolTable->entries();
    for (G4int i=0; i<nEntries; i++) {
        VolTableEntry *vte = VolTable->at(i);
        if (vte->vname == *vname) {
            return vte;
        }
    }
    return NULL;
}

G4int G3VolTable::Increment(G4String *vname)
{
    curEntry = find(vname);
    if (curEntry == NULL) {
        return 0;
    } else {
        curEntry->count++;
        return curEntry->count;
    }
}
