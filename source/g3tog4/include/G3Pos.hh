#ifndef _G3POS_
#define _G3POS_ 1

#include "G4ThreeVector.hh"
#include "G3VolTable.hh"

class G3Pos{
private:
  VolTableEntry* _VTE;
  G4int _Copy;
  G4String _LV;
  G4String _Mother;
  VolTableEntry* _VTEMother;
  G4ThreeVector* _Position;
  G4int _Irot;
  G4String _Only;
public:
  G3Pos(const G4String& LV, const G4int C, const G4String& M, 
	G4ThreeVector* T, const G4int R, const G4String& O);

  G4bool operator == (const G3Pos& g3p) const;

  virtual ~G3Pos(){;}

  VolTableEntry* GetMotherVTE();

  G4String GetName();

  G4String GetMotherName();
};
#endif












