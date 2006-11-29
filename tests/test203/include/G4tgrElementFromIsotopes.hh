#ifndef G4tgrElementFromIsotopes_h
#define G4tgrElementFromIsotopes_h

using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrElementFromIsotopes    
Author:      P. Arce
Changes:     14/07/00: creation  
----------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Transient class of a chemical element. 
*/

#include <vector>
#include "G4tgrElement.hh"

//--------------------------------------------------------------------  
class G4tgrElementFromIsotopes : public G4tgrElement { 

 public:    
  G4tgrElementFromIsotopes(){ };
  ~G4tgrElementFromIsotopes(){ };

  //! construct the G4tgrElementFromIsotopes (fill its data members) interpreting the data in the list of words 'wl' 
  G4tgrElementFromIsotopes( const vector<G4String>& wl );

  const G4int GetNumberOfIsotopes() const {return theNoIsotopes;}
  const G4String GetComponent( G4int n ) {return theComponents[n];}
  const G4double GetAbundance( G4int n ) {return theAbundances[n];}

 private:
  G4int theNoIsotopes;
  vector<G4String>  theComponents;
  vector<G4double>  theAbundances;

};

#endif

