#ifndef B01MaterialFactory_hh
#define B01MaterialFactory_hh B01MaterialFactory_hh

#include "globals.hh"
#include "g4std/map"

class G4Material;
class G4Element;

typedef G4std::map< G4String , G4Element* > B01MapSymbolElement;
typedef G4std::map< G4Element* , G4double > B01MapElementFraction;

class B01MaterialFactory{
public:
  B01MaterialFactory();
  ~B01MaterialFactory();
  
  G4Material *CreateConcrete();
  G4Material *CreateLightConcrete();
  G4Material *CreateGalactic();
  
private:
  B01MapSymbolElement fMapSymbolElement;
  B01MapElementFraction fConcreteFractions;
  void FillElementMap(const G4String &name, 
		      const G4String &symbol,
		      G4int Z,
		      G4double A);
};


#endif
