#ifndef G4tgrIsotope_h
#define G4tgrIsotope_h

using namespace std;
#include "globals.hh"
#include <vector>

/*---------------------------------------------------------------------------   
ClassName:   G4tgrIsotope    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a chemical element. 
*/  


//----------------------------------------------------------------------------  
class G4tgrIsotope { 

 public:    
  G4tgrIsotope(){ };
  ~G4tgrIsotope(){ };

  //! construct the G4tgrIsotope (fill its data members) interpreting the data in the list of words 'wl' 
  G4tgrIsotope( const vector<G4String>& wl );
  
  // Retrieval methods
  G4String GetName()  const {return theName;};
  G4int    GetZ()     const {return theZ;};
  G4int    GetN()     const {return theN;};
  G4double GetA()     const {return theA;};
  
 private:
  G4String theName;              // name of the Isotope
  G4int    theZ;                 // atomic number
  G4int    theN;                 // number of nucleons
  G4double theA;                 // mass of a 

};

#endif

