#ifndef G4tgbMaterialMixtureByWeight_h
#define G4tgbMaterialMixtureByWeight_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialMixtureByWeight    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material whose components are given by weight fractions independent of G4
 */

#include <vector>
#include <string>
#include "G4tgbMaterialMixture.hh"

class G4tgbMaterialMixtureByWeight: public G4tgbMaterialMixture {
  
public:
  G4tgbMaterialMixtureByWeight( G4tgrMaterial* tgr );
  G4tgbMaterialMixtureByWeight(){};
  ~G4tgbMaterialMixtureByWeight(){};

  //! return the associated G4Material and if does not exits, build it
  virtual G4Material* BuildG4Material();

};
#endif
