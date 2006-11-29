#ifndef G4tgbMaterialMixture_h
#define G4tgbMaterialMixture_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbMaterialMixture    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material independent of G4
 */

#include "G4tgbMaterial.hh"
#include <vector>
#include <string>

class G4tgbMaterialMixture : public G4tgbMaterial{

public:

  G4tgbMaterialMixture(){};
  ~G4tgbMaterialMixture(){};

  //Get methods
  virtual G4String GetComponent(int i) const{
    return theTgrMate->GetComponent( i );} 
  virtual double GetFraction(int i) const{
    return theTgrMate->GetFraction( i );} 

  //Operators
  //-  G4bool       operator==(const CMSMaterialMixture&) const; //Compares ONLY names
  //-  G4bool       operator!=(const CMSMaterialMixture&) const; //Compares ONLY names
  G4tgbMaterialMixture& operator= (const G4tgbMaterialMixture&); 

protected:
  virtual void TransformToFractionsByWeight(){};

};
#endif
