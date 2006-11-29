#ifndef G4tgrMaterial_h
#define G4tgrMaterial_h

using namespace std;
#include "globals.hh"
#include <iostream>
#include <vector>

/*---------------------------------------------------------------------------   
ClassName:   G4tgrMaterial    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to represent a material 
 */


class G4tgrMaterial {

  //  friend std::ostream& operator<<(std::ostream&, const G4tgrMaterial&);
  
public:

  G4tgrMaterial(){};
  virtual ~G4tgrMaterial(){};


  //Get methods
  const G4String GetName() const {return theName;}  
  //! density in g/cm3. 
  double GetDensity() const {return theDensity;}
  int GetNumberOfComponents() const {return theNoComponents;} 
  const G4String GetType() const {return theMateType;}  
  virtual double GetA() const = 0;
  virtual double GetZ() const = 0;

  virtual G4String GetComponent(int i) = 0;
  virtual double GetFraction(int i) = 0;

  G4tgrMaterial& operator= (const G4tgrMaterial&); 

protected:
  G4String  theName;
  double  theDensity;
  int  theNoComponents;
  G4String  theMateType;
};
#endif
