//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4tgbVolumeMgr.cc 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgbVolumeMgr

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbVolumeMgr.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgbVolume.hh"
#include "G4tgbMaterialMgr.hh"
#include "G4tgbRotationMatrixMgr.hh"

#include "G4tgrVolumeMgr.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrUtils.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ReflectionFactory.hh"
#include "G4tgrMessenger.hh"
#include "G4tgbDetectorBuilder.hh"

G4ThreadLocal G4tgbVolumeMgr* G4tgbVolumeMgr::theInstance = 0;

//---------------------------------------------------------------------
G4tgbVolumeMgr::G4tgbVolumeMgr() 
{
  G4ReflectionFactory::Instance()->SetScalePrecision(1.E-6*mm);
    // NOTE: problems building matrices with not enough figures,
    // like  :ROTM RR30 0.866025 0.5 0. -0.5 0.866025 0. 0. 0 -1
  theDetectorBuilder = new G4tgbDetectorBuilder();
}


//---------------------------------------------------------------------
G4tgbVolumeMgr::~G4tgbVolumeMgr()
{
  delete theDetectorBuilder;
  delete theInstance;
}


//---------------------------------------------------------------------
G4tgbVolumeMgr* G4tgbVolumeMgr::GetInstance()
{
  if( !theInstance )
  {
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
  const G4tgrVolume* tgrVoltop = theDetectorBuilder->ReadDetector();
  return theDetectorBuilder->ConstructDetector(tgrVoltop);
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4tgbVolume* vol )
{
  theVolumeList.insert( G4mssvol::value_type( vol->GetName(),
                        const_cast<G4tgbVolume*>(vol) ) );
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4VSolid* solid )
{
  theSolids.insert( G4mmssol::value_type( solid->GetName(),
                    const_cast<G4VSolid*>(solid) ) );
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4LogicalVolume* lv )
{
  theLVs.insert( G4mmslv::value_type( lv->GetName(),
                 const_cast<G4LogicalVolume*>(lv) ) );

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbVolumeMgr::RegisterMe() - Logical volume registered: "
           << lv->GetName() << G4endl; 
  }
#endif
}

//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterMe( const G4VPhysicalVolume* pv )
{
  thePVs.insert( G4mmspv::value_type( pv->GetName(),
                 const_cast<G4VPhysicalVolume*>(pv) ) );
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::RegisterChildParentLVs( const G4LogicalVolume* logvol,
                                             const G4LogicalVolume* parentLV )
{
  theLVInvTree[const_cast<G4LogicalVolume*>(logvol)] =
    const_cast<G4LogicalVolume*>(parentLV);
  theLVTree[const_cast<G4LogicalVolume*>(parentLV)] =
    const_cast<G4LogicalVolume*>(logvol);
}

//---------------------------------------------------------------------
void G4tgbVolumeMgr::CopyVolumes()
{
  //--------- Loop G4tgbVolume's and create a G4tgbVolume for each DetUnit
  G4mapsvol::iterator cite;
  G4mapsvol vollist = G4tgrVolumeMgr::GetInstance()->GetVolumeMap();
  for(cite = vollist.begin(); cite != vollist.end(); cite++)
  {
    G4tgrVolume* tgrvol = const_cast<G4tgrVolume*>( (*cite).second );
    G4tgbVolume* svol = new G4tgbVolume( tgrvol );
    RegisterMe( svol );
  }
}


//---------------------------------------------------------------------
G4tgbVolume* G4tgbVolumeMgr::FindVolume( const G4String& volname)
{
  G4mssvol::const_iterator cite = theVolumeList.find( volname );
  if( cite == theVolumeList.end() )
  {
    G4String ErrMessage = "G4tgbVolume not found: " + volname + " !";
    G4Exception("G4tgbVolumeMgr::FindVolume()", "InvalidSetup",
                FatalException, ErrMessage);
  }
  return (*cite).second;
}


//---------------------------------------------------------------------
G4VSolid* G4tgbVolumeMgr::FindG4Solid( const G4String& name )
{
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbVolumeMgr::FindG4Solid() - " << name << G4endl;
  }
#endif

  G4VSolid* oldSolid = 0;
  std::pair<G4mmssol::iterator, G4mmssol::iterator> mmssdi;
  mmssdi = theSolids.equal_range( name );

  if( mmssdi.first != mmssdi.second ) { // check there is a solid found
    G4mmssol::const_iterator mmsscite = mmssdi.first;

#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgbVolumeMgr::FindG4Solid() - Solid finding "
	     << name << G4endl; 
    }
#endif
    /*
       G4VSolid overwrites the operator== comparing the addresses
       => this cannot be used !!
       Then just compare solids by name =>> POSP tag cannot be used
       for the moment ...
         if( solid == *( (*mmsscite).second ) )
         {
           oldSolid = (*mmsscite).second;
           break;
         }
       until we write operator== for each solid type, we take a solid
       with the same name (therefore we will not allow two solids with
       equal name and different parameters (POSP) )
    */
    oldSolid = (*mmsscite).second;
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
    {
      G4cout << " G4tgbVolumeMgr::FindG4Solid() - Solid already found "
	     << name << G4endl; 
    }
#endif
  }
 
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
      G4cout << " G4tgbVolumeMgr::FindG4Solid() - Old solid: "
	     << oldSolid << G4endl;
  }
#endif

  return oldSolid;
}


//---------------------------------------------------------------------
G4LogicalVolume*
G4tgbVolumeMgr::FindG4LogVol( const G4String& name, const G4bool exists )
{
  G4mmslv::const_iterator mscite = theLVs.find( name );
  if( mscite == theLVs.end() )
  {
    if( exists )
    {
      G4String ErrMessage = "Logical volume name " + name + " not found !";
      G4Exception("G4tgbVolumeMgr::FindG4LogVol()", "InvalidSetup",
                  FatalException, ErrMessage);
    }
    return 0;
  }
  else
  {
    return (*mscite).second;
  }
}

//---------------------------------------------------------------------
G4VPhysicalVolume*
G4tgbVolumeMgr::FindG4PhysVol( const G4String& name, const G4bool exists )
{
  G4mmspv::const_iterator mscite = thePVs.find( name );
  if( mscite == thePVs.end() )
  {
    if( exists )
    {
      G4String ErrMessage = "Physical volume name " + name + " not found !";
      G4Exception("G4tgbVolumeMgr::FindG4PhysVol()", "InvalidSetup",
                  FatalException, ErrMessage);
    }
    return 0;
  }
  else
  {
    return (*mscite).second;
  }
}


//---------------------------------------------------------------------
G4VPhysicalVolume* G4tgbVolumeMgr::GetTopPhysVol()
{
  G4LogicalVolume* lv = GetTopLogVol(); 
  G4VPhysicalVolume* pv = ( *(thePVs.find( lv->GetName() )) ).second;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbVolumeMgr::GetTopPhysVol() - pv: "
           << pv->GetName() << G4endl;
  }
#endif

  return pv;
}


//---------------------------------------------------------------------
G4LogicalVolume* G4tgbVolumeMgr::GetTopLogVol()
{
  //----------- Start from any G4LogicalVolume, because if you go upwards
  //            you will always end at the top  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgbVolumeMgr::GetTopLogVol theLVInvTresize "
           << theLVInvTree.size() << G4endl;
  }
#endif
  if(  theLVInvTree.size() == 0 ) 
  { 
    G4Exception("G4tgbVolumeMgr::GetTopLogVol()", "InvalidSetup",
                FatalException, "theLVInvTree has no elements.");
  }
  G4LogicalVolume* lv = (*(theLVInvTree.begin())).second;

  //------- if first element is the top LV, its parent is 0
  if( lv == 0 ) 
  {
    lv = (*(theLVInvTree.begin())).first;
  } 
  else 
  {
    while( (*(theLVInvTree.find( lv ))).second != 0) 
    {
      //---------- get parent of first position
      lv = (*(theLVInvTree.find( lv ))).second;
#ifdef G4VERBOSE
      if( G4tgrMessenger::GetVerboseLevel() >= 2 )
      {
        G4cout << " G4tgbVolumeMgr::GetTopPhysVol: lv "
               << lv->GetName() << G4endl;
      }
#endif
    }
  }
  
  return lv;
}

  
//---------------------------------------------------------------------
void G4tgbVolumeMgr::BuildPhysVolTree()
{
/*
  G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*>::iterator cite;
  for( cite = pvs->begin(); cite != pvs->end(); cite++ )
  {
    thePVTree[ *cite ] = (*cite)->GetMother();
    thePVInvTree[ (*cite)->GetMother() ] = *cite;
  }
*/
}
 

//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpSummary()
{
  //---------- Dump number of objects of each class
  G4cout << " @@@@@@@@@@@@@ Dumping Geant4 geometry objects Summary " << G4endl;
  G4cout << " @@@ Geometry built inside world volume: "
         << GetTopPhysVol()->GetName() << G4endl;
  G4cout << " Number of G4VSolid's: " << theSolids.size() << G4endl;
  G4cout << " Number of G4LogicalVolume's: " << theLVs.size() << G4endl;
  G4cout << " Number of G4VPhysicalVolume's: " << thePVs.size() << G4endl;

  G4tgbMaterialMgr* mateMgr = G4tgbMaterialMgr::GetInstance();
  G4cout << " Number of G4Isotope's: "
         << mateMgr->GetG4IsotopeList().size() << G4endl;
  G4cout << " Number of G4Element's: "
         << mateMgr->GetG4ElementList().size() << G4endl;
  G4cout << " Number of G4Material's: "
         << mateMgr->GetG4MaterialList().size() << G4endl;

  G4tgbRotationMatrixMgr* rotmMgr = G4tgbRotationMatrixMgr::GetInstance();
  G4cout << " Number of G4RotationMatrix's: "
         << rotmMgr->GetG4RotMatList().size() << G4endl;

  //---------- Dump list of objects of each class
  DumpG4SolidList();
  DumpG4LogVolTree();
  DumpG4PhysVolTree();
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4SolidList()
{
  G4mmssol::const_iterator cite;
  for( cite = theSolids.begin(); cite != theSolids.end(); cite++)
  {
    G4cout << "G4SOLID: " << (*cite).second->GetName()
           << " of type " << (*cite).second->GetEntityType() << G4endl;
  }
} 


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4LogVolTree()
{
  G4cout << " @@@@@@@@@@@@@ DUMPING G4LogicalVolume's Tree  " << G4endl;
  
  G4LogicalVolume* lv = GetTopLogVol();
  
  DumpG4LogVolLeaf(lv, 0);
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4LogVolLeaf( const G4LogicalVolume* lv,
                                       unsigned int leafDepth)
{
  for( size_t ii=0; ii < leafDepth; ii++ )
  {
    G4cout << "  ";
  }
  G4cout << " LV:(" << leafDepth << ")" << lv->GetName() << G4endl;
  
  //---------- construct the children of this volume
  // G4LogicalVolume* lvnc = const_cast<G4LogicalVolume*>(lv);
  // std::pair<G4mlvlv::iterator, G4mlvlv::iterator> children
  //   = theLVTree.equal_range( lvnc );
  //
  // G4mlvlv::iterator cite; 
  
  leafDepth++;
  // for( cite = children.first; cite != children.second; cite++ )
  // {
  //  DumpG4LVLeaf( (*cite)->second, leafDepth );
  // } 
}



//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4PhysVolTree()
{
  G4cout << " @@@@@@@@@@@@@ DUMPING G4PhysicalVolume's Tree  " << G4endl;

  G4VPhysicalVolume* pv = GetTopPhysVol();

  DumpG4PhysVolLeaf(pv, 0);
}


//---------------------------------------------------------------------
void G4tgbVolumeMgr::DumpG4PhysVolLeaf( const G4VPhysicalVolume* pv,
                                        unsigned int leafDepth)
{
  for( size_t ii=0; ii < leafDepth; ii++ )
  {
    G4cout << "  ";
  }
  G4cout << " PV:(" << leafDepth << ")" << pv->GetName() << G4endl;

  //---------- construct the children of this PV
  // G4VPhysicalVolume* pvnc = const_cast<G4VPhysicalVolume*>(pv);
  // std::pair<G4mpvpv::iterator, G4mpvpv::iterator> children
  //  = thePVTree.equal_range( pvnc );
  //
  // G4mpvpv::iterator cite; 

  leafDepth++;
  // for( cite = children.first; cite != children.second; cite++ )
  // {
  //  DumpG4PVLeaf( (*cite)->second, leafDepth );
  // }
}
