#ifndef TiaraMaterials_hh
#define TiaraMaterials_hh TiaraMaterials_hh

#include "globals.hh"
#include "g4std/map"

class G4Material;
class G4Element;

typedef G4std::map< G4String , G4Element* > TiaraMapSymbolElement;
typedef G4std::map< G4String, G4Material* > TiaraMapNameMaterial;

class TiaraMaterials{
public:
  TiaraMaterials();
  ~TiaraMaterials();

  
  G4Material *GetMaterial(const G4String &matName) const;
  
  G4Material *CreateAir();
  G4Material *CreateConcrete();
  G4Material *CreateMCNPConcrete();
  G4Material *CreateIron();
  G4Material *CreateVakuum();

private:
  TiaraMaterials(const TiaraMaterials &);

  void FillElementMap(const G4String &name, 
		      const G4String &symbol,
		      G4int Z,
		      G4double A);

  TiaraMaterials &operator=(const TiaraMaterials &);

  TiaraMapSymbolElement fMapSymbolElement;
  TiaraMapNameMaterial fMapNameMaterial;
};

#endif
