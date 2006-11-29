#ifndef G4tgrElement_h
#define G4tgrElement_h

using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrElement    
Author:      P. Arce
Changes:     14/07/00: creation  
----------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a chemical element. 
*/

#include <vector>

//--------------------------------------------------------------------  
class G4tgrElement { 

 public:    
  G4tgrElement(){ };
  ~G4tgrElement(){ };

  const G4String& GetName() const{ return theName; }
  const G4String& GetSymbol() const{ return theSymbol; }
  const G4String& GetType() const{ return theType; }

 protected:
  G4String theName;
  G4String theSymbol;
  G4String theType;
};

#endif

