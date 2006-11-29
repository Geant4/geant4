#ifndef G4tgbElement_h
#define G4tgbElement_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbElement    
Author:      P. Arce
Changes:     14/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a chemical element. 
Build a G4Element
*/  

#include <vector>
#include <string>
#include "G4tgrElement.hh"
#include "G4Element.hh"

//----------------------------------------------------------------------------  
class G4tgbElement { 

 public:    
  G4tgbElement(){ };
  ~G4tgbElement(){ };

  //! construct the G4tgbElement from the corresponding HgElement
  G4tgbElement( G4tgrElement* hg );
  
  //! build a G4Element out of this HElement
  G4Element* BuildG4ElementSimple( );
  G4Element* BuildG4ElementFromIsotopes();

 public:
  const G4String GetName() const{
    return theTgrElem->GetName();
  }
  const G4String GetType() const{
    return theTgrElem->GetType();
  }

 private:
  G4tgrElement* theTgrElem;
  G4Element* theG4Elem;

};

#endif

