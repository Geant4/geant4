#ifndef Tst33MaterialFactory_hh
#define Tst33MaterialFactory_hh Tst33MaterialFactory_hh

#include "globals.hh"
#include "g4std/map"

class G4Material;
class G4Element;

typedef G4std::map< G4String , G4Element* > Tst33MapSymbolElement;
typedef G4std::map< G4Element* , G4double > Tst33MapElementFraction;

class Tst33MaterialFactory{
public:
  Tst33MaterialFactory();
  ~Tst33MaterialFactory();
  
  G4Material *CreateConcrete();
  G4Material *CreateLightConcrete();
  G4Material *CreateGalactic();
  
private:
  Tst33MaterialFactory(const Tst33MaterialFactory &);

  void FillElementMap(const G4String &name, 
		      const G4String &symbol,
		      G4int Z,
		      G4double A);

  Tst33MaterialFactory &operator=(const Tst33MaterialFactory &);

  Tst33MapSymbolElement fMapSymbolElement;
  Tst33MapElementFraction fConcreteFractions;
};


#endif
