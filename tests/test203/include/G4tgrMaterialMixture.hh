#ifndef G4tgrMaterialMixture_h
#define G4tgrMaterialMixture_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrMaterialMixture    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material mixture
 */

#include "G4tgrMaterial.hh"
#include <vector>

class G4tgrMaterialMixture : public G4tgrMaterial{

  friend ostream& operator<<(ostream&, const G4tgrMaterialMixture&);
  
public:

  G4tgrMaterialMixture(){};
  ~G4tgrMaterialMixture(){};

  //! fill the data interpreting the list of words read 'wl'
  G4tgrMaterialMixture(const G4String& matType, const vector<G4String>& wl);

  //Get methods
  virtual double GetA() const {return 0.;}
  virtual double GetZ() const {return 0.;}
  virtual G4String GetComponent(int i) {return theComponents[i];} 
  virtual double GetFraction(int i) {return theFractions[i];} 

  G4tgrMaterialMixture operator= (const G4tgrMaterialMixture&); 

protected:
  virtual void TransformToFractionsByWeight(){};

protected:
  vector<G4String>  theComponents;
  vector<double>  theFractions;
};
#endif
