#ifndef _G3VOLTABLE_
#define _G3VOLTABLE_ 1

#include <rw/tphdict.h>

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VSolid;

class VolTableEntry {
private:
  G4String _Vname;
  G4String _Shape;
  G4double* _Rpar;
  G4int _Npar;
  G4Material* _Mat;
  G4VSolid* _Solid;
  G4LogicalVolume* _LV;
  G4int _Code;
  friend class G3VolTable;

public:
  VolTableEntry(G4String& v, G4String& sh, G4double* R, G4int n, G4Material* m,
		G4VSolid* so, G4LogicalVolume* lv, G4int code);
  virtual ~VolTableEntry();
  G4bool NegVolPars();
  G4bool Deferred();
};
  
class G3VolTable{
private:
  G4LogicalVolume* G3toG4LogicalMother;
  G4VPhysicalVolume* _pv;
  G4int _copy;
  G4VPhysicalVolume* G3toG4PhysicalMother;
  RWTPtrHashDictionary <G4String, VolTableEntry>* _VTD;
  VolTableEntry* _VTE;

private:
  void SetMother(G4LogicalVolume* theMother);

public:
  void PutLV(G4String& vname, G4String& shape, G4double* Rpar, G4int Npar,
	     G4Material* mat, G4VSolid* solid, G4LogicalVolume* lv, G4int cd);

  VolTableEntry* GetVTE(const G4String& Vname);

  G4LogicalVolume* GetLV(const G4String& Vname, G4String& Shape, 
			 G4double*& Rpar, G4int& Npar, G4Material*& Mat, 
			 G4VSolid*& Solid, G4int& Code);

  G4LogicalVolume* GetLV(const G4String& vname="-");

  G4VPhysicalVolume* GetPV(const G4String& pname="-", G4int pcopy=0);

  G3VolTable();

  virtual ~G3VolTable();
};
extern G3VolTable G3Vol;
#endif












