#ifndef _VOLTABLEENTRY_
#define _VOLTABLEENTRY_ 1

#include <rw/tphdict.h>
#include <rw/tpordvec.h>
#include "G4ThreeVector.hh"
#include "G3Pos.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class VolTableEntry {
private:
  G4int _Nmed;
  G4String _Vname;
  G4String _Shape;
  G4double* _Rpar;
  G4int _Npar;
  G4Material* _Mat;
  G4VSolid* _Solid;
  G4LogicalVolume* _LV;
  G4bool _Deferred;
  G4bool _NegVolPars;
  RWTPtrOrderedVector <VolTableEntry> _Daughters; // VolTableEntry Daughters
  RWTPtrOrderedVector <G3Pos> _G3Pos; // associated G3Pos objects

public:
  VolTableEntry(G4String& v, G4String& sh, G4double* R, 
		G4int n, G4int nmed, G4Material* m, 
		G4VSolid* so, G4bool Deferred, 
		G4bool NegVolPars);

  virtual ~VolTableEntry();

  void SetLV(G4LogicalVolume* ll);

  G4String GetName();

  G4VSolid* GetSolid();

  void AddG3Pos(G3Pos* aG3Pos);

  G4int NPCopies();

  G3Pos* GetG3PosCopy(G4int copy=0);

  G4bool HasNegVolPars();

  G4bool HasDeferred();

  G4String GetShape();

  G4double* GetRpar();

  G4int GetNpar();

  G4Material* GetMaterial();

  G4LogicalVolume* GetLV();

  G4int GetNmed();

  void AddDaughter(VolTableEntry* aDaughter);

  G4int GetNoDaughters();

  VolTableEntry* GetDaughter(G4int i);

  VolTableEntry* FindDaughter(const G4String& vname);

  G4bool operator == ( const VolTableEntry& vte) const;
};
#endif

