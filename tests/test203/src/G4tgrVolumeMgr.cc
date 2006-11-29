// G4tgrVolumeMgr
#include "G4tgrVolumeMgr.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMaterialFactory.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrSolid.hh"
#include "G4tgrSolidBoolean.hh"


G4tgrVolumeMgr* G4tgrVolumeMgr::theInstance = 0;

//-------------------------------------------------------------
G4tgrVolumeMgr* G4tgrVolumeMgr::GetInstance()
{
  if( !theInstance ) {
    theInstance = new G4tgrVolumeMgr;
  }
  return theInstance;
}


/*//-------------------------------------------------------------
void G4tgrVolumeMgr::fillDU( const vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  if( wl.size() < 5 ) {
    G4tgrUtils::dumpVS(wl, "!!!! EXITING: G4tgrVolumeMgr::addDU. Line read with less than 5 words ");
    exit(1);
  }

  G4tgrVolume* du = findDU( G4tgrUtils::subQuotes( wl[1] ) );
  if( du == 0 ) {
    new G4tgrVolume();
    du->constructVolume( wl );
  } else {
    cerr << "!!!! EXITING: G4tgrVolume already exists: tag :VOLU repeated " << G4tgrUtils::subQuotes( wl[1] ) << G4endl;
    exit(1);
  }

}*/



//-------------------------------------------------------------------
G4tgrSolid* G4tgrVolumeMgr::CreateSolid( const std::vector<G4String>& wl, G4bool bVOLUtag )
{
  G4tgrSolid* sol = FindSolid( wl[1] );
  if( sol ) {
    G4Exception("G4tgrVolumeMgr::CreateSolid  Solid already exists " + wl[1]);
  } 
  
  std::vector<G4String> wlc = wl;
  if( bVOLUtag ) wlc.pop_back();

  //-  G4cout << " G4tgrVolumeMgr::CreateSolid isFromVOLU " << bVOLUtag << " " << wlc.size() << G4endl;
  
  G4String wl2 = wlc[2];
  for( size_t ii = 0; ii < wl2.length(); ii++ ){
    wl2[ii] = toupper( wl2[ii] );
  }
  if( wl2 == "UNION" || wl2 == "SUBS" || wl2 == "INTERS" ) {
    //-------------------------------- boolean solid
    //---------- Create G4tgrSolidBoolean and fill the solid params
    sol = new G4tgrSolidBoolean( wlc );
  } else {
    //---------- Create G4tgrSolidSimple and fill the solid params
    sol = new G4tgrSolid( wlc );
  }

  //  G4cout << " G4tgrVolumeMgr::CreateSolid " << sol << " " << sol->GetName() << G4endl;

  return sol;
}

//-------------------------------------------------------------------
void G4tgrVolumeMgr::RegisterMe( G4tgrSolid* sol) 
{
  if( theG4tgrSolidMap.find( sol->GetName() ) != theG4tgrSolidMap.end() ) { 
    G4Exception(" G4tgrVolumeMgr: cannot be two solids with the same name " + sol->GetName());
  }
  theG4tgrSolidMap.insert(mapssol::value_type(sol->GetName(), sol) ); 
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::UnRegisterMe( G4tgrSolid* sol ) 
{
  if( theG4tgrSolidMap.find( sol->GetName() ) != theG4tgrSolidMap.end() ) { 
    G4Exception("G4tgrSolidMgr::unRegisterMe: cannot unregister a Solid that is not registered " + sol->GetName() );
  } else {
    theG4tgrSolidMap.erase( theG4tgrSolidMap.find( sol->GetName() ) ); 
  }
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::RegisterMe( G4tgrVolume* vol) 
{
  theG4tgrVolumeList.push_back( vol );
  if( theG4tgrVolumeMap.find( vol->GetName() ) != theG4tgrVolumeMap.end() ) { 
    G4Exception(" G4tgrVolumeMgr: cannot be two volumes with the same name " + vol->GetName());
  }
  theG4tgrVolumeMap.insert(mapsvol::value_type(vol->GetName(), vol) ); 
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::UnRegisterMe( G4tgrVolume* vol ) 
{
  vector<G4tgrVolume*>::iterator ite;
  for(ite = theG4tgrVolumeList.begin(); ite != theG4tgrVolumeList.end(); ite++) {
    if((*ite) == vol ) break;
  }
  if( ite ==  theG4tgrVolumeList.end() ) { 
    G4Exception("G4tgrVolumeMgr::unRegisterMe: cannot unregister a Volume that is not registered " + vol->GetName() );
  } else {
    theG4tgrVolumeList.erase( ite );
  }
  theG4tgrVolumeMap.erase( theG4tgrVolumeMap.find( vol->GetName() ) ); 
}

//-------------------------------------------------------------

void G4tgrVolumeMgr::RegisterParentChild( const G4String& parentName, const G4tgrPlace* pl )
{ 
  theG4tgrVolumeTree.insert(mmapspl::value_type(parentName, pl) );

}


//-------------------------------------------------------------
G4tgrSolid* G4tgrVolumeMgr::FindSolid( const G4String& volname, bool exists ) 
{
  G4tgrSolid* vol = 0;
  
  mapssol::iterator svite = theG4tgrSolidMap.find( volname );
  if( svite == theG4tgrSolidMap.end() ) {
    if( exists ) {
      for( svite = theG4tgrSolidMap.begin();  svite != theG4tgrSolidMap.end(); svite++ ) {
	cerr << " VOL:" << (*svite).first << G4endl;
      }
      G4Exception("!!!! EXITING: G4tgrSolid not found in list: " + volname );
    }
  } else {
    vol = const_cast<G4tgrSolid*>((*svite).second);
  }

  return vol;
}

//-------------------------------------------------------------
G4tgrVolume* G4tgrVolumeMgr::FindVolume( const G4String& volname, bool exists ) 
{
  G4tgrVolume* vol = 0;
  
  mapsvol::iterator svite = theG4tgrVolumeMap.find( volname );
  if( svite == theG4tgrVolumeMap.end() ) {
    if( exists ) {
      for( svite = theG4tgrVolumeMap.begin();  svite != theG4tgrVolumeMap.end(); svite++ ) {
	cerr << " VOL:" << (*svite).first << G4endl;
      }
      G4Exception("!!!! EXITING: G4tgrVolume not found in list: " + volname );
    }
  } else {
    vol = const_cast<G4tgrVolume*>((*svite).second);
  }

  return vol;
}


//-------------------------------------------------------------
const G4tgrVolume* G4tgrVolumeMgr::GetTopVolume()
{
  //----------- Start from any G4tgrVolume, because if you go upwards you will always end at the top  
  const G4tgrVolume* vol = (*(theG4tgrVolumeMap.begin())).second;
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 ) 
    G4cout << " G4tgrVolumeMgr::GetTopVolume: vol " << vol->GetName() << " no place = " <<  vol->GetPlacements().size() << G4endl;
#endif

  //-  G4String vol1stParent = *(vol->placements().begin())->parentName();
  while( vol->GetPlacements().size() != 0 ) {
    //---------- get parent of first placement (there could be a pathological case and you will
    //---------- never get to the top, but it wll not be the case for HARP
    vol = FindVolume( (* (vol->GetPlacements()).begin() )->GetParentName() );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << " G4tgrVolumeMgr::GetTopVolume: vol " << vol->GetName()<< " no place = " <<  vol->GetPlacements().size() << G4endl;
#endif
  };

  return vol;
}


//-------------------------------------------------------------
pair<mmapspl::iterator, mmapspl::iterator> G4tgrVolumeMgr::GetChildren( const G4String& name )
{
  pair<mmapspl::iterator, mmapspl::iterator> dite;
  dite = theG4tgrVolumeTree.equal_range( name );
  return dite;
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::DumpVolumeTree()
{
  G4cout << " @@@@@@@@@@@@@@@@ DUMPING G4tgrVolume's Tree  " << G4endl;

  const G4tgrVolume* vol = GetTopVolume();

  DumpVolumeLeaf( vol, 0, 0);
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::DumpVolumeLeaf( const G4tgrVolume* vol, uint copyNo, uint leafDepth)
{
  for( uint ii=0; ii < leafDepth; ii++ ){
    G4cout << "  ";
  }
  G4cout << " VOL:(" << leafDepth << ")" <<  vol->GetName() << "   copy No " << copyNo << G4endl;

  //---------- construct the children of this VOL
  pair<mmapspl::const_iterator, mmapspl::const_iterator> children = 
   GetChildren( vol->GetName() );
  mmapspl::const_iterator cite; 

  leafDepth++;
  for( cite = children.first; cite != children.second; cite++ ) {
    //---- find G4tgrVolume pointed by G4tgrPlace
    const G4tgrPlace* pla = (*cite).second;
    const G4tgrVolume* volchild = pla->GetVolume();
    //--- find copyNo
    uint cn = pla->GetCopyNo();
    DumpVolumeLeaf( volchild, cn, leafDepth );
  }
}


//-------------------------------------------------------------
void G4tgrVolumeMgr::DumpSummary()
{
  //---------- Dump number of objects of each class
  G4cout << " @@@@@@@@@@@@@@@@@@ Dumping Detector Summary " << G4endl;
  G4cout << " @@@ Geometry build inside world volume: " << GetTopVolume()->GetName() << G4endl;
  G4cout << " Number of G4tgrVolume's: " << theG4tgrVolumeMap.size() << G4endl;
  mapsvol::const_iterator cite;
  uint nPlace = 0;
  for( cite = theG4tgrVolumeMap.begin();cite != theG4tgrVolumeMap.end(); cite++ ) {
    nPlace += ((*cite).second)->GetPlacements().size();
  }
  G4cout << " Number of G4tgrPlace's: " << nPlace << G4endl;

  G4tgrMaterialFactory* matef = G4tgrMaterialFactory::GetInstance();
  G4cout << " Number of G4tgrIsotope's: " << matef->GetIsotopeList().size() << G4endl;
  G4cout << " Number of G4tgrElement's: " << matef->GetElementList().size() << G4endl;
  G4cout << " Number of G4tgrMaterial's: " << matef->GetMaterialList().size() << G4endl;

  G4tgrRotationMatrixFactory* rotmf = G4tgrRotationMatrixFactory::GetInstance();
  G4cout << " Number of G4tgrRotationMatrix's: " << rotmf->GetRotMatList().size() << G4endl;


  //---------- Dump detail list of objects of each class
  DumpVolumeTree();
  
  matef->DumpIsotopeList();
  matef->DumpElementList();
  matef->DumpMaterialList();

  rotmf->DumpRotmList();
}
