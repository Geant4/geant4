// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3VolTableEntry.hh,v 1.3 1999-12-09 01:27:46 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

class G3VolTableEntry {
  public:
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











