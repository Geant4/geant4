#include "G4tgrVolume.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlaceSimple.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"
#include <set>

#include "CLHEP/Units/SystemOfUnits.h"


//-------------------------------------------------------------
G4tgrVolume::G4tgrVolume( const vector<G4String>& wl) 
{
  theType = "VOLSimple";

  //---------- set name 
  theName = G4tgrUtils::SubQuotes( wl[1] ); 

  if( wl.size() != 4 ) {
    //:VOLU tag to build a volume creating solid and material
    //---------- set material name
    theMaterialName = G4tgrUtils::SubQuotes( wl[wl.size()-1] );
    
    //---------- create only vector<double> of theSolidParams
    theSolid = FindOrCreateSolid( wl );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << "G4tgrVolume::G4tgrVolume  from new  solid " << theSolid->GetName() << " of type " << theSolid->GetType() << G4endl;
#endif
  } else {
    //:VOLU tag to build a volume assigning material to solid
    //---------- set material name
    theMaterialName = G4tgrUtils::SubQuotes( wl[3] );

    theSolid = G4tgrVolumeMgr::GetInstance()->FindSolid( wl[2], TRUE );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << "G4tgrVolume::G4tgrVolume  from existing solid " << theSolid->GetName() << " of type " << theSolid->GetType() << G4endl;
#endif
  }

  theVisibility = 1;
  theRGBColour = new double[3];
  for(size_t ii=0; ii<3; ii++) theRGBColour[ii] = -1.;
}


//-------------------------------------------------------------------------
G4tgrSolid* G4tgrVolume::FindOrCreateSolid( const std::vector<G4String>& wl )
{
  G4tgrSolid* solid = G4tgrVolumeMgr::GetInstance()->FindSolid( wl[2] );
  if( !solid ) {
    solid = G4tgrVolumeMgr::GetInstance()->CreateSolid( wl, 1 );
  }
  //  G4cout << " G4tgrVolume::CreateSolid " << solid->GetName() << G4endl;
  return solid;

}

//-------------------------------------------------------------
G4tgrPlace* G4tgrVolume::AddPlace( const vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 8, WLSIZE_EQ, " G4tgrVolume::AddPlace");
  //---------- set G4tgrPlace 
  G4tgrPlaceSimple* pl = new G4tgrPlaceSimple( wl );
  //---------- check that there is no previous placement in the same parent with the same copyNo 
  vector<G4tgrPlace*>::iterator ite;
  for( ite = thePlacements.begin(); ite != thePlacements.end(); ite++) {
    if( (*ite)->GetCopyNo() == pl->GetCopyNo() && (*ite)->GetParentName() == pl->GetParentName()  ) { 
      cerr << "G4tgrVolume::AddPlace repeated placement " << pl->GetCopyNo() << " for volume" << GetName() << G4endl;
      exit(1);
    }
  } 

  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << pl << " new placement " << thePlacements.size() << " added for Volume " << theName << " inside " << pl->GetParentName() << " type " << pl->GetType() << G4endl;
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}



//-------------------------------------------------------------
G4tgrPlaceDivRep* G4tgrVolume::AddPlaceReplica( const std::vector<G4String>& wl ) 
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE, " G4tgrVolume::AddPlaceReplica");
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_LE, " G4tgrVolume::AddPlaceReplica");

  //---------- set G4tgrPlace 
  G4tgrPlaceDivRep* pl = new G4tgrPlaceDivRep( wl );
  pl->SetType("PlaceReplica");
  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << " new placement " << thePlacements.size() << " added for Volume " << theName << " inside " << pl->GetParentName() << G4endl;
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}


//-------------------------------------------------------------
G4tgrPlaceParameterisation* G4tgrVolume::AddPlaceParam( const vector<G4String>& wl )
{
  //---------- set G4tgrPlaceParameterisation
  G4tgrPlaceParameterisation* pl = new G4tgrPlaceParameterisation( wl );

  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    G4cout << " new placement Param " << thePlacements.size() << " added for Volume " << theName << " inside " << pl->GetParentName() << G4endl;
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}


//-------------------------------------------------------------
void G4tgrVolume::AddVisibility( const vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 3, WLSIZE_EQ, " G4tgrVolume::AddVisibility");

  //---------- set G4tgrPlace 
  theVisibility = G4tgrUtils::GetBool( wl[2] ); 

}


//-------------------------------------------------------------
void G4tgrVolume::AddRGBColour( const vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ, " G4tgrVolume::AddRGBColour");
  
  //  theRGBColour = new double[3];
  theRGBColour[0] = G4tgrUtils::GetFloat( wl[2] );
  theRGBColour[1] = G4tgrUtils::GetFloat( wl[3] );
  theRGBColour[2] = G4tgrUtils::GetFloat( wl[4] );
}



