#ifndef G4tgrVolumeMgr_h
#define G4tgrVolumeMgr_h
#include "globals.hh"

/*---------------------------------------------------------------------------   
ClassName:   G4tgrVolumeMgr    
Author:      P. Arce
Changes:     03/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the detector units.
It is a singleton, accesed always with GetInstance() */ 

#include "G4tgrSolid.hh"
#include "G4tgrVolume.hh"
#include "G4tgrPlace.hh"
#include "G4tgrIsotope.hh"
#include "G4tgrElement.hh"
#include "G4tgrMaterial.hh"
#include "G4tgrRotationMatrix.hh"

#include <map>
//#include <mmap>
//#include <algo>

typedef map< const G4String, const G4tgrSolid* > mapssol;
typedef map< const G4String, const G4tgrVolume* > mapsvol;
typedef multimap< const G4String, const G4tgrPlace* > mmapspl;

//----------------------------------------------------------------------------  
class G4tgrVolumeMgr 
{ 
 public:    

  G4tgrVolumeMgr(){ };
  ~G4tgrVolumeMgr(){ 
    //delete theG4tgrVolumeList;
    //delete theG4tgrVolumeTree;
  };


  /// Get the only instance 
  static G4tgrVolumeMgr* GetInstance();  

  G4tgrSolid* CreateSolid( const std::vector<G4String>& wl, G4bool bVOLUtag ); 

  //! add to theG4tgrVolumeTree
  void RegisterParentChild( const G4String& parentName, const G4tgrPlace* pl );

  //! find a G4tgrSolid with name 'name'. If it is not found: if exists is true , exits; if exists is false, return 0
  G4tgrSolid* FindSolid( const G4String& name, bool exists = false );

  //! find a G4tgrVolume with name 'volname'. If it is not found: if exists is true , exits; if exists is false, return 0
  G4tgrVolume* FindVolume( const G4String& volname, bool exists = false );

  //! find the top of the volume tree
  const G4tgrVolume* GetTopVolume();  

  //! find the list of G4tgrPlace children of  G4tgrVolume 'name'
  pair<mmapspl::iterator, mmapspl::iterator> GetChildren( const G4String& name );

  //! dump summary
  void DumpSummary();
  //! dump to cout the tree of G4tgrVolume's
  void DumpVolumeTree();
  //! dump a G4tgrVolume indicating its copy no and its depth (=number of antecesors)
  void DumpVolumeLeaf( const G4tgrVolume* vol, uint copyNo, uint leafDepth);

  //! for GAUDI
  void RegisterMe( G4tgrSolid* vol);
  void UnRegisterMe( G4tgrSolid* vol );
  void RegisterMe( G4tgrVolume* vol);
  void UnRegisterMe( G4tgrVolume* vol );
  void RegisterMe( G4tgrPlace* pl ) { theG4tgrPlaceList.push_back( pl ); }
  void RegisterMe( G4tgrIsotope* iso ) { theHgIsotList.push_back( iso ); }
  void RegisterMe( G4tgrElement* ele ) { theHgElemList.push_back( ele ); }
  void RegisterMe( G4tgrMaterial* mat ) { theHgMateList.push_back( mat ); }
  void RegisterMe( G4tgrRotationMatrix* rm ) { theHgRotMList.push_back( rm ); }

//! public data access functions
 public:
  const mapssol& GetSolidMap() {return theG4tgrSolidMap;}
  const mapsvol& GetVolumeMap() {return theG4tgrVolumeMap;}
  const mmapspl& GetVolumeTree() {return theG4tgrVolumeTree;}
  vector<G4tgrVolume*> GetVolumeList() {return theG4tgrVolumeList;}
  vector<G4tgrPlace*> GetDetPlaceList() {return theG4tgrPlaceList;}
  vector<G4tgrIsotope*> GetIsotopeList() {return theHgIsotList;}
  vector<G4tgrElement*> GetElementList() {return theHgElemList;}
  vector<G4tgrMaterial*> GetMaterialList() {return theHgMateList;}
  vector<G4tgrRotationMatrix*> GetRotMList() {return theHgRotMList;}

private:
  //! map of G4tgrSolid's: G4String is the G4tgrSolid name, G4tgrSolid* the pointer to it
  mapssol theG4tgrSolidMap;
  //! map of G4tgrVolume's: G4String is the G4tgrVolume name, G4tgrVolume* the pointer to it
  mapsvol theG4tgrVolumeMap;

  //! hierarchy tree of G4tgrVolume's: G4String is the name of the parent G4tgrVolume, G4tgrPlace* the pointers to children
  mmapspl theG4tgrVolumeTree;

  static G4tgrVolumeMgr* theInstance;

  vector<G4tgrVolume*> theG4tgrVolumeList;
  vector<G4tgrPlace*> theG4tgrPlaceList;
  vector<G4tgrIsotope*> theHgIsotList;
  vector<G4tgrElement*> theHgElemList;
  vector<G4tgrMaterial*> theHgMateList;
  vector<G4tgrRotationMatrix*> theHgRotMList;

};

#endif
