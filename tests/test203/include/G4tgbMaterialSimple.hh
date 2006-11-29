#ifndef G4tgbMaterialSimple_h
#define G4tgbMaterialSimple_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialSimple    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material independent of G4
 */

#include <vector>
#include <string>
#include "G4tgbMaterial.hh"


class G4tgbMaterialSimple : public G4tgbMaterial {

 public:
 
  G4tgbMaterialSimple(){};
  ~G4tgbMaterialSimple(){};

  //! fill the data interpreting the list of words read 'wl'
  G4tgbMaterialSimple( G4tgrMaterial* tgr );

  friend ostream& operator<<(ostream&, const G4tgbMaterialSimple&);

  //! return the associated G4Material and if does not exits, build it
  virtual G4Material* BuildG4Material();
  G4double GetZ() const {return theZ;}
  G4double GetA() const {return theA;}

private:
  G4double    theZ;
  G4double    theA;

};
#endif
