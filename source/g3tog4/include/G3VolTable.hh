#ifndef _G3VOLTABLE_
#define _G3VOLTABLE_ 1

#include <rw/tphdict.h>
#include "VolTableEntry.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class G3VolTable{
private:
  VolTableEntry* G3toG4TopVTE;
  G4String _FirstKey;
  VolTableEntry* _VTE;
  RWTPtrHashDictionary <G4String, VolTableEntry>* _VTD;
  G4int _NVTE;
  G4int _NG3Pos;

public:
  void PutVTE(G4String& vname, G4String& shape, 
	      G4double* Rpar,  G4int Npar, G4int nmed,
	      G4Material* mat, G4VSolid* solid, 
	      G4bool Deferred, G4bool ng);
  
  RWTPtrHashDictionary <G4String, VolTableEntry>* GetVTD() ;
  VolTableEntry* GetVTE(G4String& Vname);
  G3VolTable();
  virtual ~G3VolTable();
  G4LogicalVolume* GetG3toG4Mother();
  VolTableEntry* GetFirstVTE();
  void SetFirstVTE();
  void VTEStat();
  void CountG3Pos();
};
extern G3VolTable G3Vol;
#endif












