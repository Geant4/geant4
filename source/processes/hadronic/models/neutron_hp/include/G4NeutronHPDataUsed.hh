#ifndef G4NeutronHPDataUsed_h
#define G4NeutronHPDataUsed_h 1

#include "globals.hh"

class G4NeutronHPDataUsed
{
  public:
  
  G4NeutronHPDataUsed();
  
  void SetA(G4double anA){theA = anA;}
  void SetZ(G4int aZ){theZ = aZ;}
  void SetName(G4String aName){theName = aName;}

  G4int GetZ() {return theZ;}
  G4double GetA() {return theA;}
  G4String GetName() {return theName;}
  
  private:
  
  G4String theName;
  G4double theA;
  G4int theZ;
};

#endif
