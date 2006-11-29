#include "G4tgbVolume.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgbRotationMatrixMgr.hh"

#include "G4tgrVolumeMgr.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrUtils.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ReflectionFactory.hh"
#include "G4tgrMessenger.hh"

G4tgbVolumeMgr* G4tgbVolumeMgr::theInstance = 0;

//---------------------------------------------------------------------
G4tgbVolumeMgr::G4tgbVolumeMgr() 
{
  G4ReflectionFactory::Instance()->SetScalePrecision(1.E-6*mm); // problems building matrices with not enough figures, like  :ROTM RR30 0.866025 0.5 0. -0.5 0.866025 0. 0. 0 -1
}


//---------------------------------------------------------------------
G4tgbVolumeMgr* G4tgbVolumeMgr::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgbVolumeMgr();
  }
  return theInstance;
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::AddTextFile( const G4String& fname ) 
{
  G4tgrFileReader::GetInstance()->AddTextFile( fname );
}


//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbVolumeMgr::ReadAndConstructDetector()
{
  const G4tgrVolume* tgrVoltop = ReadDetector();
  return ConstructDetector(tgrVoltop);
}


//---------------------------------------------------------------------
const G4tgrVolume* G4tgbVolumeMgr::ReadDetector()
{
  //------------------- construct g4 geometry
  //---------- find top G4tgrVolume 
  G4tgrFileReader::GetInstance()->Initialize();
  G4tgrVolumeMgr* tgrVolmgr = G4tgrVolumeMgr::GetInstance();
  const G4tgrVolume* tgrVoltop = tgrVolmgr->GetTopVolume();  
  return tgrVoltop;
}


//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbVolumeMgr::ConstructDetector( const G4tgrVolume* tgrVoltop )
{

  //---------- copy list of G4tgrVolume's to list of G4tgbVolume's (just a trick to make all GEANT4 volume building in this class)
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  tgbVolmgr->CopyVolumes();
  //---------- find corresponding volume in list of G4tgbVolume's
  G4tgbVolume* tgbVoltop = tgbVolmgr->FindVolume( tgrVoltop->GetName() );
  
  //---------- ConstructG4Volumes of top G4tgbVolume (it will recursively build the whole tree)
  tgbVoltop->ConstructG4Volumes( 0, (const G4LogicalVolume*)0 );

  G4VPhysicalVolume* physvol = (G4tgbVolumeMgr::GetInstance())->GetTopPhysVol();
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
  cout << " G4tgbDetectorConstruction::Construct. TopPV " << physvol->GetName() << endl;
#endif
  return physvol;

}

//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4tgbVolume* vol )
{
  theVolumeList.insert( mssvol::value_type( vol->GetName(), const_cast<G4tgbVolume*>(vol) ) );

}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4VSolid* solid )
{
  //  cout << " registering sg4vsolid " << solid->GetName() << endl;
  theSolids.insert( mmssol::value_type( solid->GetName(), const_cast<G4VSolid*>(solid) ) );

}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4LogicalVolume* lv )
{
  theLVs.insert( mmslv::value_type( lv->GetName(), const_cast<G4LogicalVolume*>(lv) ) );

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::RegisterMe. log vol registered " << lv->GetName() << endl; 
#endif
}

//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4VPhysicalVolume* pv )
{
  thePVs.insert( mmspv::value_type( pv->GetName(), const_cast<G4VPhysicalVolume*>(pv) ) );

}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::CopyVolumes()
{
  //--------- Loop G4tgbVolume's and create a G4tgbVolume for each DetUnit
  mapsvol::iterator cite;
  mapsvol vollist = G4tgrVolumeMgr::GetInstance()->GetVolumeMap();
  for(cite = vollist.begin(); cite != vollist.end(); cite++) {
    G4tgrVolume* tgrvol = const_cast<G4tgrVolume*>( (*cite).second );
    G4tgbVolume* svol = new G4tgbVolume( tgrvol );
    RegisterMe( svol );
  }
 
}


//---------------------------------------------------------------------
G4tgbVolume* G4tgbVolumeMgr::FindVolume( const G4String& volname)
{
  mssvol::const_iterator cite = theVolumeList.find( volname );
  if( cite == theVolumeList.end() ) {
    G4Exception("!!!! EXITING: G4tgbVolume not found: " + volname);
  }
  return (*cite).second;

}


//---------------------------------------------------------------------
G4VSolid* G4tgbVolumeMgr::FindG4Solid( const G4String& name )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << "  G4tgbVolumeMgr::FindG4Solid " << name << endl;
#endif

  G4VSolid* oldSolid = 0;
  mmssol::const_iterator mmsscite;
  pair<mmssol::iterator, mmssol::iterator> mmssdi;
  mmssdi = theSolids.equal_range( name );

  //-  cout << " solids equal range size " << &(*(mmssdi.first)) << "  "  << &(*(mmssdi.second)) << " solids.size " << theSolids.size() << endl;
  for( mmsscite = mmssdi.first; mmsscite != mmssdi.second; mmsscite++ ) {
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::FindG4Solid: solid finding " << name << endl; 
#endif
    /* !! G4VSolid overwrites the operator== comparing the addresses => this cannot be used !!
Then just compare solids by name =>> POSP tag cannot be used for the moment
    if( solid == *( (*mmsscite).second ) ) {
      oldSolid = (*mmsscite).second;
      if( verbose >= 2) cout << "  G4tgbVolumeMgr::FindG4Solid: solid already found " << solid.GetName()  << endl; 
      break;
      } */ 
    // until we write operator== for each solid type, we take a solid with the same name (therefore we will not allow two solids with equal name and different parameters (POSP) )
    oldSolid = (*mmsscite).second;
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
    cout << "  G4tgbVolumeMgr::FindG4Solid: solid already found " << name << endl; 
#endif
    break;
  }
 
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::FindG4Solid oldsolid: " << oldSolid << endl;
#endif
  return oldSolid;
}


//---------------------------------------------------------------------
G4LogicalVolume* G4tgbVolumeMgr::FindG4LogVol( const G4String& name, const G4bool exists )
{
  //-cout << " G4tgbVolumeMgr::FindLV " << name << endl;

  mmslv::const_iterator mscite = theLVs.find( name );
 //- cout << " lvs equal range size " << &(*(mmsdi.first)) << "  "  << &(*(mmsdi.second)) << " LVs.size " << theLVs.size() << endl;
  if( mscite == theLVs.end() ) {
    if( exists ) {
      G4Exception("!!! EXITING:  G4tgbVolumeMgr::FindLV. Logical Volume name " + name + " not found ");
    }
    return 0;
  } else {
    return (*mscite).second;
  }

}


//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbVolumeMgr::GetTopPhysVol()
{
  G4LogicalVolume* lv = GetTopLogVol(); 
  G4VPhysicalVolume* pv = ( *(thePVs.find( lv->GetName() )) ).second;
  //  cout << "lv " << lv->GetName() << " pv " << pv << endl;
  // const G4VPhysicalVolume* pv = ( *(thePVs.find( "expHall" )) ).second;
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::GetTopPhysVol: pv " << pv->GetName() << endl;
#endif

  return pv;
  
}


//---------------------------------------------------------------------
G4LogicalVolume* G4tgbVolumeMgr::GetTopLogVol()
{
  //----------- Start from any G4LogicalVolume, because if you go upwards you will always end at the top  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::GetTopLogVol theLVInvTresize " << theLVInvTree.size() << endl;
#endif
  if(  theLVInvTree.size() == 0 ) { 
    G4Exception("!!! EXITING:  G4tgbVolumeMgr::GetTopLogVol. theLVInvTree has no elements");
  }
  G4LogicalVolume* lv = (*(theLVInvTree.begin())).second;
  //------- if first element is the top LV, its parent is 0
  if( lv == 0 ) {
    lv = (*(theLVInvTree.begin())).first;
  } else {
    while( (*(theLVInvTree.find( lv ))).second != 0) {
      //---------- get parent of first position
      lv = (*(theLVInvTree.find( lv ))).second;
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 ) 
    cout << " G4tgbVolumeMgr::GetTopPhysVol: lv " << lv->GetName() << endl;
#endif
    };
  }

  return lv;
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::BuildPhysVolTree()
{
//  if( G4RunManager::GetGeometryManager()->is initialised )

  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  vector<G4VPhysicalVolume*>::iterator cite;
  for( cite = pvs->begin(); cite != pvs->end(); cite++ ) {
    //t    thePVTree[ *cite ] = (*cite)->GetMother();
    //t     thePVInvTree[ (*cite)->GetMother() ] = *cite;
  }

}
 

//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpSummary()
{
  //---------- Dump number of objects of each class
  cout << " @@@@@@@@@@@@@@@@@@ Dumping G4 geometry objects Summary " << endl;
  cout << " @@@ Geometry build inside world volume: " << GetTopPhysVol()->GetName() << endl;
  cout << " Number of G4VSolid's: " << theSolids.size() << endl;
  cout << " Number of G4LogicalVolume's: " << theLVs.size() << endl;
  cout << " Number of G4VPhysicalVolume's: " << thePVs.size() << endl;

  G4tgbMaterialMgr* mateMgr = G4tgbMaterialMgr::GetInstance();
  cout << " Number of G4Isotope's: " << mateMgr->GetG4IsotopeList().size() << endl;
  cout << " Number of G4Element's: " << mateMgr->GetG4ElementList().size() << endl;
  cout << " Number of G4Material's: " << mateMgr->GetG4MaterialList().size() << endl;

  G4tgbRotationMatrixMgr* rotmMgr = G4tgbRotationMatrixMgr::GetInstance();
  cout << " Number of G4RotationMatrix's: " << rotmMgr->GetG4RotMatList().size() << endl;

 
  //---------- Dump list of objects of each class
  DumpG4SolidList();
  DumpG4LogVolTree();
  DumpG4PhysVolTree();
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4SolidList()
{
  mmssol::const_iterator cite;
  for( cite = theSolids.begin(); cite != theSolids.end(); cite++) {
    cout << "G4SOLID: " << (*cite).second->GetName() << " of type " << (*cite).second->GetEntityType() << endl;
  }
} 


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4LogVolTree()
{
  cout << " @@@@@@@@@@@@@@@@ DUMPING G4LogicalVolume's Tree  " << endl;
  
  G4LogicalVolume* lv = GetTopLogVol();
  
  DumpG4LogVolLeaf( lv, 0);
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4LogVolLeaf( const G4LogicalVolume* lv, uint leafDepth)
{
  for( uint ii=0; ii < leafDepth; ii++ ){
    cout << "  ";
  }
  cout << " LV:(" << leafDepth << ")" << lv->GetName() << endl;
  
  //---------- construct the children of this volume
  G4LogicalVolume* lvnc = const_cast<G4LogicalVolume*>(lv);
  pair<mlvlv::iterator, mlvlv::iterator> children = theLVTree.equal_range( lvnc );
  
  mlvlv::iterator cite; 
  
  leafDepth++;
   for( cite = children.first; cite != children.second; cite++ ) {
     //t       DumpG4LVLeaf( (*cite)->second, leafDepth );
    } 
}



//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4PhysVolTree()
{
  cout << " @@@@@@@@@@@@@@@@ DUMPING G4PhysicalVolume's Tree  " << endl;

  G4VPhysicalVolume* pv = GetTopPhysVol();

  DumpG4PhysVolLeaf( pv, 0);
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4PhysVolLeaf( const G4VPhysicalVolume* pv, uint leafDepth)
{
  for( uint ii=0; ii < leafDepth; ii++ ){
    cout << "  ";
  }
  cout << " PV:(" << leafDepth << ")" << pv->GetName() << endl;

  //---------- construct the children of this PV
  G4VPhysicalVolume* pvnc = const_cast<G4VPhysicalVolume*>(pv);
  pair<mpvpv::iterator, mpvpv::iterator> children = thePVTree.equal_range( pvnc );

  mpvpv::iterator cite; 

  leafDepth++;
  for( cite = children.first; cite != children.second; cite++ ) {
    //t    DumpG4PVLeaf( (*cite)->second, leafDepth );
  }
}
