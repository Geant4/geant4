#ifndef G4tgbMaterial_h
#define G4tgbMaterial_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterial    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material independent of G4
 */

#include <vector>
#include <string>
#include "G4tgrMaterial.hh"

#include "G4Material.hh"

class G4tgbMaterial {

  friend ostream& operator<<(ostream&, const G4tgbMaterial&);
  
 public:

  G4tgbMaterial(){};
  virtual ~G4tgbMaterial(){};
  G4tgbMaterial( G4tgrMaterial* tgr );

  virtual G4Material* BuildG4Material() = 0;

  //!
 public:
  const G4String GetName() const {
    return theTgrMate->GetName();
  }
  const double GetDensity() const {
    return theTgrMate->GetDensity();
  }
  const int GetNumberOfMaterials() const {
    return theTgrMate->GetNumberOfComponents();
  }
  const double GetA() const {
    return theTgrMate->GetA();
  }
  const double GetZ() const {
    return theTgrMate->GetZ();
  }

  const G4String GetType() const {
    return theTgrMate->GetType();
  }  

 protected:
  G4tgrMaterial* theTgrMate;
  G4Material* theG4Mate;

};
#endif
