#ifndef G4tgrPlace_h
#define G4tgrPlace_h

using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrPlace
Author:      P. Arce
Changes:     03/07/00: creation  
             01/01: make it abstract base class
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Abstract base class to descripe the position of a G4tgrVolume inside another G4tgrVolume */ 
//----------------------------------------------------------------------------  

#include <vector>
#include "CLHEP/Vector/ThreeVector.h"
class G4tgrVolume;


class G4tgrPlace
{
 public:
  G4tgrPlace(){ };
  virtual ~G4tgrPlace(){ };

  //acces private data members
 public:
  const G4String GetParentName() const {return theParentName;}

  G4tgrVolume* GetVolume() const {return theVolume;}

  uint GetCopyNo() const {return theCopyNo;}

  G4String GetType() const {return theType;}

  void SetVolume( G4tgrVolume* vol ) { theVolume = vol; }

  void SetType( const G4String& typ ){ theType = typ; }

  virtual G4String GetRotMatName() const { return G4String(" "); };
  virtual CLHEP::Hep3Vector GetPlacement() const { return (CLHEP::Hep3Vector)0; }; //??

 protected:
  //! the detunit to which it belongs
  G4tgrVolume* theVolume;

  //! the parent (by name, as we will allow that a child is placed in the file before the parent is created) 
  G4String theParentName;  

  //! the copy number
  uint theCopyNo;

  //! the type (simple/replica/param)
  G4String theType;

};
#endif

