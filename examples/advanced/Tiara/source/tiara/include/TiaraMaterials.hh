// $Id: TiaraMaterials.hh,v 1.3 2003-06-18 16:40:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraMaterials
//

#ifndef TiaraMaterials_hh
#define TiaraMaterials_hh TiaraMaterials_hh

#include "globals.hh"
#include <map>

class G4Material;
class G4Element;

typedef std::map< G4String , G4Element* > TiaraMapSymbolElement;
typedef std::map< G4String, G4Material* > TiaraMapNameMaterial;

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
