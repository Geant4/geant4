#ifndef _G3VOLTABLE_
#define _G3VOLTABLE_ 1

#include "g4rw/tphdict.h"
#include "VolTableEntry.hh"

class G4LogicalVolume;
class G4Material;
class G4VSolid;

class G3VolTable{
private:
  VolTableEntry* G3toG4TopVTE;
  G4String _FirstKey;
  VolTableEntry* _VTE;
  G4RWTPtrHashDictionary <G4String, VolTableEntry>* _VTD;
  G4int _NVTE;
  G4int _NG3Pos;

public:
  VolTableEntry* PutVTE(VolTableEntry* aVTE);  
  G4RWTPtrHashDictionary <G4String, VolTableEntry>* GetVTD() ;
  VolTableEntry* GetVTE(const G4String& Vname);
  void ListVTE();
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












