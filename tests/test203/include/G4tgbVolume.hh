#ifndef G4tgbVolume_h
#define G4tgbVolume_h
using namespace std;
#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgbVolume    
Author:      P. Arce
Changes:     03/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the geometry info of any detector unit (independent of G4) */ 
/*! The detector units created in this class are essentially transient copies of G4 
physical volumes. Thus, they are characterized by a name and the parameters of a 
G4 physical volume. 

They have associated several detector positions, that can be instances of G4tgrPlace, G4tgrPlaceDivRep or G4tgrPlaceParameterisation  
Each detector positioning is done inside a parent. As there can be several parents, 
we will write one parent for each det. pos., even if that means that parents 
are repeated

*/  

#include <vector>
#include <map>
#include "G4tgrVolume.hh"
class G4tgrPlace;
class G4tgrSolid;

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

//----------------------------------------------------------------------------  
class G4tgbVolume { 

 public:    
  G4tgbVolume(){};
  ~G4tgbVolume(){ };
  G4tgbVolume( G4tgrVolume* vol);

  //! construct the G4VSolid, G4LogicalVolume and the G4VPhysicalVolume of copy 'copyNo'
  void ConstructG4Volumes( const G4tgrPlace* place, const G4LogicalVolume* parentLV );

  //! construct the G4VSolid from the data of the corresponding G4tgrVolume. 
  //! Allow to use data from another G4tgrVolume, needed by boolean solids (that have to construct two solids and then do the boolean operation)
  G4VSolid* FindOrConstructG4Solid( const G4tgrSolid* vol);

  //! construct the G4LogicalVolume and then call the construction of volumes that are positioned inside this LV
  G4LogicalVolume* ConstructG4LogVol( const G4VSolid* solid );

  //! construct the G4VPhysicalVolume placing curentLV with position given by the G4tgrPlace 'copyNo' inside 'parentLV'
  G4VPhysicalVolume* ConstructG4PhysVol( const G4tgrPlace* place, const G4LogicalVolume* currentLV, const G4LogicalVolume* parentLV  );


  void SetCutsInRange( G4LogicalVolume* logvol, std::map<G4String,double> cuts );
  void SetCutsInEnergy( G4LogicalVolume* logvol, std::map<G4String,double> cuts );


  //! before building a solid of type 'solydType', check if the number of paramenters is the expected one
  void CheckNoSolidParams( const G4String& solidType, const int NoParamExpected, const uint NoParam );

 public:
  const G4String GetName() const{
    return theTgrVolume->GetName();
  }
  const G4bool GetVisibility() const {
    return theTgrVolume->GetVisibility();
  }
  const G4double* GetColour() const {
    return theTgrVolume->GetColour();
  }

//! private data members
 private:   
  //! the G4tgrVolume to which it corresponds
  G4tgrVolume* theTgrVolume;
};

#endif
 
