#ifndef G4tgrElementSimple_h
#define G4tgrElementSimple_h

using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrElementSimple    
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
class G4tgrElementSimple : public G4tgrElement{ 

 public:    
  G4tgrElementSimple(){ };
  ~G4tgrElementSimple(){ };

  //! construct the G4tgrElementSimple (fill its data members) interpreting the data in the list of words 'wl' 
  G4tgrElementSimple( const vector<G4String>& wl );

  const double GetZ() const{ return theZ; }
  const double GetA() const{ return theA; }

 private:
  double theZ;
  double theA;

};

#endif

