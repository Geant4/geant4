#ifndef _G3POS_
#define _G3POS_ 1

#include "G4ThreeVector.hh"

class VolTableEntry;

class G3Pos{
private:
  VolTableEntry* _VTE;
  G4int _Copy;
  G4String _LV;
  G4ThreeVector* _Position;
  G4int _Irot;
  G4String _Only;
public:
  G3Pos(){;}

  G3Pos(G4String& LV, G4int C, G4ThreeVector* T, G4int R, G4String& O);

  G4bool operator == (const G3Pos& g3p) const;

  virtual ~G3Pos();

  VolTableEntry* GetVTE();

  G4String& GetName();

  G4int GetIrot();

  G4ThreeVector* GetPos();

  G4int GetCopy();
};
#endif

