#ifndef G4tgrVolume_h
#define G4tgrVolume_h

#include "globals.hh"
/*---------------------------------------------------------------------------   
ClassName:   G4tgrVolume    
Author:      P. Arce
Changes:     03/07/00: creation  
             01/01: make it abstract base class
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Abstract base class to manage the geometry info of any volume  */ 
/*! The volumes created in this class contain the information of a detector 
volume.
They have associated several detector placements, that can be instances of G4tgrPlace, G4tgrPlaceDivision, G4tgrPlaceDivRep or G4tgrPlaceParameterisation  
Each detector positioning is done inside a parent. As there can be several parents, 
we will write one parent for each volume placement, even if that means that parents 
are repeated
*/  

#include <vector>
#include <map>
#include "globals.hh"
#include "G4tgrSolid.hh"
#include "G4tgrPlace.hh"
class G4tgrPlaceDivRep;
class G4tgrPlaceParameterisation;
//----------------------------------------------------------------------------  
class G4tgrVolume { 

 public:    

  G4tgrVolume(){};
  G4tgrVolume( const vector<G4String>& wl );
  virtual ~G4tgrVolume(){ };

  // Find a solid and if it does no exist, create it
  G4tgrSolid* FindOrCreateSolid( const std::vector<G4String>& wl );

  //! add a position with the data read from a ':pos' tag
  G4tgrPlace* AddPlace( const vector<G4String>& wl );

  //! add a replicated position
  G4tgrPlaceDivRep* AddPlaceReplica( const vector<G4String>& wl );

  //! add a parameterised position
  G4tgrPlaceParameterisation* AddPlaceParam( const vector<G4String>& wl );

  //! add visibility flag
  void AddVisibility( const vector<G4String>& wl );

  //! add colour
  void AddRGBColour( const vector<G4String>& wl );

  /*
  //! add a replica position (a division) with the data read from a ':dvn'or a ':dvn2' tag
  bool buildDivisionByNoDiv( vector<G4String> wl );

  //! add a replica position (a division) with the data read from a ':dvt'or a ':dvt2' tag
  bool buildDivisionByDivWidth( vector<G4String> wl );
  */

  //! Public methods to get and set private data members  
  const G4String GetName() const {return theName;}
  const G4String GetType() const {return theType;}
  const G4tgrSolid* GetSolid() const {return theSolid;}
  const G4String GetMaterialName() const {return theMaterialName;}

  const vector<G4tgrPlace*> GetPlacements() const {return thePlacements;}
  const G4bool GetVisibility() const {return theVisibility;}
  const G4double* GetColour() const {return theRGBColour;}
  virtual G4tgrVolume* GetVolume( int ii ) const{
    cerr << " !!EXITING: G4tgrVolume::getDU. should only be called for G4tgrVolumeWithBooleanSolid objects " << ii << endl;
    exit(1);
    return (G4tgrVolume*)(0);
  }

//! private data members
 protected:   
  //! Name of the volume
  G4String theName;   
  //! Type of the volume    
  G4String theType;   
  //! Material of which the corresponding PV will be made of.   
  G4String theMaterialName;   
  //! Solid 
  G4tgrSolid* theSolid;

 //! Vector of placements 
  vector<G4tgrPlace*> thePlacements;
    
 protected:
  bool theVisibility;
 
  double* theRGBColour;

};

#endif

