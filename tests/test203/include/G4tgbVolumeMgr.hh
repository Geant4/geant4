#ifndef G4tgbVolumeMgr_h
#define G4tgbVolumeMgr_h
#include "globals.hh"

/*---------------------------------------------------------------------------   
Author:      P. Arce
Changes:     03/08/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class to manage the G4 volumes: G4VSolids, G4LogicalVolumes, G4VPhysicalVolumes
It is a singleton, accesed always with GetInstance() */ 

#include <string>
#include <vector>
#include <map>

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
using namespace std;

class G4tgbVolume;
class G4tgrVolume;
 
typedef map< G4String, G4tgbVolume* > mssvol;
typedef multimap< G4String, G4VSolid* > mmssol;
typedef multimap< G4String, G4LogicalVolume* > mmslv;
typedef multimap< G4String, G4VPhysicalVolume* > mmspv;
typedef map< G4LogicalVolume*, G4LogicalVolume* > mlvlv;
typedef map< G4VPhysicalVolume*, G4VPhysicalVolume* > mpvpv;

//----------------------------------------------------------------------------  
class G4tgbVolumeMgr 
{ 
 public:    

  G4tgbVolumeMgr();
  ~G4tgbVolumeMgr(){ };

  void AddTextFile( const G4String& fname );
  G4VPhysicalVolume* ReadAndConstructDetector();
  const G4tgrVolume* ReadDetector();
  G4VPhysicalVolume* ConstructDetector( const G4tgrVolume* tgrVoltop );

  //! build a G4tgbVolume per each G4tgbVolume
  void CopyVolumes();

  //! find a G4tgbVolume by name
  G4tgbVolume* FindVolume( const G4String& volname);

  //! register a G4tgbVolume
  void RegisterMe( const G4tgbVolume* vol );
  //! register a G4VSolid
  void RegisterMe( const G4VSolid* solid );
  //! register a G4LogicalVolume
  void RegisterMe( const G4LogicalVolume* lv );
  //! register a G4VPhysicalVolume
  void RegisterMe( const G4VPhysicalVolume* pv );
  //! register a child and its parent LV
  void RegisterChildParentLVs( const G4LogicalVolume* logvol, const G4LogicalVolume* parentLV ){
    //-    cout << " registerChildParentLVs " << logvol->GetName() << " " << parentLV << endl;
    theLVInvTree[const_cast<G4LogicalVolume*>(logvol)] = const_cast<G4LogicalVolume*>(parentLV);
    theLVTree[const_cast<G4LogicalVolume*>(parentLV)] = const_cast<G4LogicalVolume*>(logvol);
  }

  /*! find if solid already exists, comparing the name and all parameters 
  (could be checked before creating it, but it would be quite complicated,   
  because it would have to compare the parameters, and they depend on the type of solid) */
  G4VSolid* FindG4Solid( const G4String& name );
  
  /*! find if LV already exists, comparing the name and the solid 
    (with POSP it may happen that the same volume has two different solidParameter and therefore two solids will be created */
  G4LogicalVolume* FindG4LogVol( const G4String& theName, const bool bExists = 0 );
  std::vector<G4LogicalVolume*> FindLVs( const G4String& theName, const G4bool exists = 0 );

  //! get the top PV in the hierarchy tree: calls topLV, because physicalvolumes are not placed until geometry is initialized
  G4VPhysicalVolume* GetTopPhysVol();
  //! get the top LV in the hierarchy tree
  G4LogicalVolume* GetTopLogVol();


  void BuildPhysVolTree();

//volmping methods
  void DumpSummary();
  void DumpG4LogVolTree();
  void DumpG4LogVolLeaf(const G4LogicalVolume* lv, uint leafDepth);
  void DumpG4PhysVolTree();
  void DumpG4PhysVolLeaf(const G4VPhysicalVolume* pv, uint leafDepth);
  void DumpG4SolidList();

  //! public access functions
 public:
  //! Get the only instance 
  static G4tgbVolumeMgr* GetInstance();  

  const multimap< G4String, G4VSolid* >& GetSolids() const {
    return theSolids;
  }

private:
  static G4tgbVolumeMgr* theInstance;

  //! map of G4tgbVolume's: G4String is the G4tgbVolume name, G4tgbVolume* the pointer to it
  mssvol theVolumeList;

  // solids container
  mmssol theSolids;
  // logical volume container
  mmslv theLVs;
  // physical volume container
  mmspv thePVs;

  // logical volume tree for navigation (from parent to children): first is parent, then child
  mlvlv theLVTree;
  // logical volume tree for inverse navigation (from children to parent): first is child, then parent
  mlvlv theLVInvTree;

  // physical volume tree for navigation (from parent to children): first is parent, then child
  mpvpv thePVTree;
  // physical volume tree for inverse navigation (from children to parents): first is child, then parent
  mpvpv thePVInvTree;

};

#endif

