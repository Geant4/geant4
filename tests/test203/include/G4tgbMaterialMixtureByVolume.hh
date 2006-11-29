#ifndef G4tgbMaterialMixtureByVolume_h
#define G4tgbMaterialMixtureByVolume_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialMixtureByVolume    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material whose components are given by volume
 */

#include <vector>
#include <string>
#include "G4tgbMaterialMixture.hh"

class G4tgbMaterialMixtureByVolume: public G4tgbMaterialMixture {
  
public:
  G4tgbMaterialMixtureByVolume( G4tgrMaterial* tgr );
  G4tgbMaterialMixtureByVolume(){};
  ~G4tgbMaterialMixtureByVolume(){};

   //! return the associated G4Material and if does not exits, build it
  virtual G4Material* BuildG4Material();

  virtual void TransformToFractionsByWeight();

private:
  double* theFractionsByWeight;
};
#endif
