#ifndef _VOLTABLEENTRY_
#define _VOLTABLEENTRY_ 1

#include <rw/tphdict.h>
#include <rw/tpordvec.h>
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;
class G3Pos;

class VolTableEntry {
private:
  G4int _Nmed;
  G4String _Vname;
  G4String _Shape;
  G4double* _Rpar;
  G4int _Npar;
  const G4Material* _Mat;
  const G4VSolid* _Solid;
  G4LogicalVolume* _LV;
  G4bool _Deferred;
  G4bool _NegVolPars;
  RWTPtrOrderedVector <G3Pos> _Daughters; // G3Pos Daughters
  RWTPtrOrderedVector <G3Pos> _G3Pos; // associated G3Pos objects

public:
  VolTableEntry(const G4String& v, const G4String& sh, const G4double* R, 
		const G4int n, const G4int nmed, const G4Material* m, 
		const G4VSolid* so, const G4bool Deferred, 
		const G4bool NegVolPars);

  virtual ~VolTableEntry();

  void SetLV(G4LogicalVolume* ll);

  G4String GetName();

  void AddG3Pos(G3Pos* aG3Pos);

  G4int NPCopies() const;

  G3Pos* GetG3PosCopy(G4int copy=0) const;

  G4bool HasNegVolPars();

  G4bool HasDeferred();

  G4String GetShape();

  G4double* GetRpar();

  G4int GetNpar();

  const G4Material* GetMaterial();

  G4LogicalVolume* GetLV();

  G4int GetNmed();

  void AddDaughter(G3Pos* aDaughter);

  G4int GetNoDaughters() const;

  G3Pos* GetDaughter(const G4int i) const;

  G4bool operator == ( const VolTableEntry& vte) const;

};
#endif












