#ifndef G4tgbIsotope_h
#define G4tgbIsotope_h
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbIsotope    
Author:      P. Arce
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of an isotope
Build a G4Isotope
*/  

#include <vector>
#include <string>
#include "G4tgrIsotope.hh"
#include "G4Isotope.hh"

//----------------------------------------------------------------------------  
class G4tgbIsotope { 

 public:    
  G4tgbIsotope(){ };
  ~G4tgbIsotope(){ };

  //! construct the G4tgbIsotope from the corresponding HgIsotope
  G4tgbIsotope( G4tgrIsotope* hg );
  
  //! build a G4Isotope out of this HIsotope
  G4Isotope* BuildG4Isotope( );

 public:
  const G4String GetName() const{
    return theTgrIsot->GetName();
  }

 private:
  G4tgrIsotope* theTgrIsot;
  G4Isotope* theG4Isot;

};

#endif

