#ifndef _G3VOLTABLE_
#define _G3VOLTABLE_ 1

#include <rw/tphdict.h>
#include "VolTableEntry.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class G3VolTable{
private:
  G4LogicalVolume* G3toG4LogicalMother;
  G4String* _FirstKey;
  VolTableEntry* _VTE;
  RWTPtrHashDictionary <G4String, VolTableEntry>* _VTD;

public:
  void PutVTE(const G4String& vname, const G4String& shape, 
	      const G4double* Rpar,  const G4int Npar, const G4int nmed,
	      const G4Material* mat, const G4VSolid* solid, 
	      const G4bool Deferred, const G4bool ng);
  
  RWTPtrHashDictionary <G4String, VolTableEntry>* GetVTD();
  VolTableEntry* GetVTE(const G4String& Vname);
  G3VolTable();
  virtual ~G3VolTable();
  G4LogicalVolume* GetG3toG4Mother();
  VolTableEntry* GetFirstVTE();
};
extern G3VolTable G3Vol;
#endif












