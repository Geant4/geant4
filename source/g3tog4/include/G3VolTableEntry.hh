// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTableEntry.hh,v 1.4 2000-11-24 09:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// This class is a temporary representation of G3 volume.
// Its methods enables successive updating of its instances
// during the phase of filling the G3 tables (defining G3 geometry,
// eg. by parsing the G3 input via clparse.cc).
// See G3VolTable class description, too.
// 
// Data members:
//  fVname       volume name; 
//  fShape       G3 shape name;
//  fRpar        array of G3 volumes parameters;
//  fNpar        number of G3 volumes parameters;
//  fNmed        G3 tracking medium number;
//  fSolid       the G4VSolid of this volume;
//  fLV          the G4LogicalVolume;
//  fHasNegPars  true if G3 volume was defined with incomplete
//               parameters; 
//  fDaughters   vector of daughter VTEs (VTEs of volumes placed inside
//               this volume); 
//  fMothers     vector of mother VTEs (VTEs of volumes inside which this
//               volume is placed);
//  fClones      vector of clone VTEs (see explanation below);
//  fG3Pos       vector of G3 positions (G3Pos objects)
//  fDivision    G3Division object created in case the G4 volume
//               was defined as division; 
//
// Clone volumes:
// In case a G3 volume (e.g. XYZ) has changed its solid parameters
// with its new position (placement with GSPOSP) a new G3VolTableEntry
// (associated with this new solid) with a new name (XYZ_N)
// is created and registered as a clone volume in the fClones vector
// data member of its master G3VolTableEntry object. 

// ----------------------
//
// by I.Hrivnacova, 13.10.99

#ifndef G3VOLTABLEENTRY_HH
#define G3VOLTABLEENTRY_HH 1

#include "globals.hh"
#include "G3Pos.hh"
#include "G3Division.hh"
#include "g4rw/tpordvec.h"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class G3VolTableEntry
{
  public:  // with description

    G3VolTableEntry(G4String& vname, G4String& shape, G4double* rpar, 
                     G4int npar, G4int nmed, G4VSolid* solid, 
		     G4bool hasNegPars);
    virtual ~G3VolTableEntry();

    // operators
    G4bool operator == ( const G3VolTableEntry& vte) const;

    // methods
    void AddG3Pos(G3Pos* aG3Pos);
    void AddDaughter(G3VolTableEntry* aDaughter);
    void AddMother(G3VolTableEntry* aDaughter);
    void AddClone(G3VolTableEntry* aDaughter);
    void ReplaceDaughter(G3VolTableEntry* vteOld, G3VolTableEntry* vteNew);
    void ReplaceMother(G3VolTableEntry* vteOld, G3VolTableEntry* vteNew);
    G3VolTableEntry* FindDaughter(const G4String& vname);
    G3VolTableEntry* FindMother(const G4String& vname);
    G3VolTableEntry* FindClone(const G4String& vname);
    void PrintSolidInfo();

    // set methods
    void SetName(G4String name);
    void SetLV(G4LogicalVolume* lv);
    void SetSolid(G4VSolid* solid);
    void SetNmed(G4int nmed);
    void SetNRpar(G4int npar, G4double* Rpar);
    void SetDivision(G3Division* division);
    void SetHasNegPars(G4bool hasNegPars);
    void ClearG3PosCopy(G4int copy);
    void ClearDivision();
 
    // get methods
    G4String  GetName();
    G4String  GetShape();
    G4int GetNmed();
    G4int GetNpar();
    G4double* GetRpar();
    G4int NPCopies();
    G3Pos* GetG3PosCopy(G4int copy=0);
    G3Division* GetDivision();
    G4bool HasNegPars();
    G4VSolid* GetSolid();
    G4LogicalVolume* GetLV();
    G4int GetNoDaughters();
    G4int GetNoMothers();
    G4int GetNoClones();
    G3VolTableEntry* GetDaughter(G4int i);
    G3VolTableEntry* GetMother(G4int i);
    G3VolTableEntry* GetMother();  
      // return the first mother - to be removed
    G3VolTableEntry* GetClone(G4int i);
    G3VolTableEntry* GetMasterClone();

  private:

    G4String fVname;
    G4String fShape;
    G4double* fRpar;
    G4int fNpar;
    G4int fNmed;
    G4VSolid* fSolid;
    G4LogicalVolume* fLV;
    G4bool fHasNegPars;
    G4RWTPtrOrderedVector<G3VolTableEntry> fDaughters;
    G4RWTPtrOrderedVector<G3VolTableEntry> fMothers;
    G4RWTPtrOrderedVector<G3VolTableEntry> fClones;
    G4RWTPtrOrderedVector<G3Pos> fG3Pos;
    G3Division*  fDivision;
};

// inline methods

inline void G3VolTableEntry::SetDivision(G3Division* division)
{ fDivision = division; }

inline G3Division* G3VolTableEntry::GetDivision()
{ return fDivision; }

#endif
