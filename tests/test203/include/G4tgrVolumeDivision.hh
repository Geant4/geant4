#ifndef G4tgrVolumeDivision_h
#define G4tgrVolumeDivision_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrVolumeDivision
Author:      P. Arce
Changes:     12/12/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
//! Class to manage the geometry info of detector unit made by dividing a mother detector unit

//#include <multimap>
#include "G4tgrVolume.hh"
#include "G4tgrPlaceDivRep.hh"

typedef multimap< G4String, G4String > mmss;

//---------------------------------------------------------------------------- 
class G4tgrVolumeDivision : public G4tgrVolume 
{ 

 public:    
  G4tgrVolumeDivision( const vector<G4String>& wl );
  virtual ~G4tgrVolumeDivision(){ };

  //! construct a volume with the data read from a ':DIV_...' tag
  //!:DIV_NUM  "volu_name" "parent_name" "material" number_divisions "replica_axis_name" offset
  //!:DIV_STEP  "volu_name" "parent_name" "material" step "replica_axis_name" offset

  //! set the solid type and parameters dividing the mother volune
  bool SetSolid(G4tgrVolume* parentDU, bool byStep, EAxis axis, double div_step, double offset);

  //! set the list of supported axis for each solid types
  static void SetSupportedAxis();

  //! Access functions
  G4tgrPlaceDivRep* GetPlaceDivision() {
    return thePlaceDiv;
  }

//! private data members
 private:
  G4tgrPlaceDivRep* thePlaceDiv;
  static mmss theSupportedAxis; 

};

#endif
 
