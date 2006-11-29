#ifndef G4tgbMaterialMixtureByNoAtoms_h
#define G4tgbMaterialMixtureByNoAtoms_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialMixtureByNoAtoms    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material whose components are given by number of atoms independent of G4
 */

#include <vector>
#include <string>
#include "G4tgbMaterialMixture.hh"

class G4tgbMaterialMixtureByNoAtoms: public G4tgbMaterialMixture {
  
public:
  G4tgbMaterialMixtureByNoAtoms( G4tgrMaterial* tgr );
  G4tgbMaterialMixtureByNoAtoms(){};
  ~G4tgbMaterialMixtureByNoAtoms(){};

  //! return the associated G4Material and if does not exits, build it
  virtual G4Material* BuildG4Material();

  virtual void TransformToFractionsByWeight();

 private:
  bool compAreElements;
  bool compAreMaterials;
};
#endif
