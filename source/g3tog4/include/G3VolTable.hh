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
  G4bool _Deferred;
  G4bool _NegVolPars;

public:
  VolTableEntry(G4String& v, G4String& sh, G4double* R, 
		G4int n, 
		G4Material* m, G4VSolid* so, 
		G4LogicalVolume* LV, G4bool Deferred, 
		G4bool NegVolPars);
  virtual ~VolTableEntry();
  G4bool HasNegVolPars(){return _NegVolPars;}
  G4bool HasDeferred(){return _Deferred;};
  G4String GetShape(){return _Shape;}
  G4double* GetRpar(){return _Rpar;}
  G4int GetNpar(){return _Npar;}
  G4Material* GetMaterial(){return _Mat;}
  G4LogicalVolume* GetLV(){return _LV;}
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
	     G4Material* mat, G4VSolid* solid, G4bool Deferred, G4bool ng);
  VolTableEntry* GetVTE(const G4String& Vname);
  G4LogicalVolume* GetLV(const G4String& vname="-");
  G4VPhysicalVolume* GetPV(const G4String& pname="-", G4int pcopy=0);
  G3VolTable();
  virtual ~G3VolTable();
};
extern G3VolTable G3Vol;
#endif












