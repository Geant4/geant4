#ifndef G4tgrMaterialSimple_h
#define G4tgrMaterialSimple_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrMaterialSimple    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a simple material, i.e. made of a single element
 */

#include <vector>
#include "G4tgrMaterial.hh"


class G4tgrMaterialSimple : public G4tgrMaterial {

 public:
 
  G4tgrMaterialSimple(){};
  ~G4tgrMaterialSimple(){};

  //! fill the data interpreting the list of words read 'wl'
  G4tgrMaterialSimple(const G4String& matType, const vector<G4String>& wl);

  friend ostream& operator<<(ostream&, const G4tgrMaterialSimple&);

  G4tgrMaterialSimple& operator= (const G4tgrMaterialSimple&); 

 public:
  virtual double GetA() const {return theA;}
  virtual double GetZ() const {return theZ;}

  virtual G4String GetComponent(int i) { 
    cerr << "!!!EXITING: component(int i) should never be called for a MaterialSimple " << i << endl; 
    exit(1); 
  }
  virtual double GetFraction(int i) {
    cerr << "!!!EXITING: fraction(int i) should never be called for a MaterialSimple " << i << endl; 
    exit(1); 
  }


 private:
  double    theA;
  double    theZ;
};
#endif
