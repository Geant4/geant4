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
// $Id: G4tgrVolume.cc,v 1.1 2008-10-23 14:43:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4tgrVolume

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrVolume.hh"
#include "G4tgrUtils.hh"
#include "G4tgrSolid.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlace.hh"
#include "G4tgrPlaceSimple.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"

//-------------------------------------------------------------
G4tgrVolume::G4tgrVolume()
{
}


//-------------------------------------------------------------
G4tgrVolume::~G4tgrVolume()
{
}


//-------------------------------------------------------------
G4tgrVolume::G4tgrVolume( const std::vector<G4String>& wl) 
{
  theType = "VOLSimple";

  //---------- set name 
  theName = G4tgrUtils::GetString( wl[1] ); 

  if( wl.size() != 4 )
  {
    //:VOLU tag to build a volume creating solid and material
    //---------- set material name
    theMaterialName = G4tgrUtils::GetString( wl[wl.size()-1] );
    
    //---------- create only vector<double> of theSolidParams
    theSolid = FindOrCreateSolid( wl );
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << "G4tgrVolume::G4tgrVolume() - From new solid: "
           << *theSolid << G4endl;
  }
#endif
  }
  else
  {
    //:VOLU tag to build a volume assigning material to solid
    //---------- set material name
    theMaterialName = G4tgrUtils::GetString( wl[3] );
    theSolid = G4tgrVolumeMgr::GetInstance()->FindSolid( wl[2], true );

#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
      G4cout << " G4tgrVolume::G4tgrVolume() - From existing solid: "
             << *theSolid << G4endl;
    }
#endif
  }
  theVisibility = 1;
  theRGBColour = new G4double[4];
  for(size_t ii=0; ii<4; ii++)  { theRGBColour[ii] = -1.; }
}


//-------------------------------------------------------------------------
G4tgrVolume* G4tgrVolume::GetVolume( G4int ii ) const
{
  G4String ErrMessage = "Should only be called for composite solids... " + ii;
  G4Exception("G4tgrVolume::GetVolume()", "InvalidCall",
              FatalException, ErrMessage);
  return 0;
}


//-------------------------------------------------------------------------
G4tgrSolid* G4tgrVolume::FindOrCreateSolid( const std::vector<G4String>& wl )
{
  G4tgrSolid* solid = G4tgrVolumeMgr::GetInstance()->FindSolid( wl[2] );
  if( !solid )
  {
    solid = G4tgrVolumeMgr::GetInstance()->CreateSolid( wl, 1 );
  }
  return solid;
}


//-------------------------------------------------------------
G4tgrPlace* G4tgrVolume::AddPlace( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 8, WLSIZE_EQ, " G4tgrVolume::AddPlace");
  //---------- set G4tgrPlace 
  G4tgrPlaceSimple* pl = new G4tgrPlaceSimple( wl );
  //---------- check that there is no previous placement in
  //           the same parent with the same copyNo 
  std::vector<G4tgrPlace*>::iterator ite;
  for( ite = thePlacements.begin(); ite != thePlacements.end(); ite++)
  {
    if( ((*ite)->GetCopyNo() == pl->GetCopyNo())
     && ((*ite)->GetParentName() == pl->GetParentName()) )
    {
      G4String ErrMessage = "Repeated placement. Volume "
                          + theName + " in " + pl->GetParentName();
      G4Exception("G4tgrVolume::AddPlace()", "InvalidArgument",
                  FatalErrorInArgument, ErrMessage);
    }
  } 

  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << "   New placement: " << thePlacements.size()
           << " added for Volume " << theName
           << " inside " << pl->GetParentName()
           << " type " << pl->GetType() << G4endl;
  }
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}


//-------------------------------------------------------------
G4tgrPlaceDivRep*
G4tgrVolume::AddPlaceReplica( const std::vector<G4String>& wl ) 
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
  {
    G4cout << "   New placement: " << thePlacements.size()
           << " added for Volume " << theName
           << " inside " << pl->GetParentName() << G4endl;
  }
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}


//-------------------------------------------------------------
G4tgrPlaceParameterisation*
G4tgrVolume::AddPlaceParam( const std::vector<G4String>& wl )
{
  //---------- set G4tgrPlaceParameterisation
  G4tgrPlaceParameterisation* pl = new G4tgrPlaceParameterisation( wl );

  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << "   New placement Param: " << thePlacements.size()
           << " added for Volume " << theName
           << " inside " << pl->GetParentName() << G4endl;
  }
#endif
  //---------- register parent - child 
  G4tgrVolumeMgr::GetInstance()->RegisterParentChild( pl->GetParentName(), pl );

  return pl;
}


//-------------------------------------------------------------
void G4tgrVolume::AddVisibility( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 3, WLSIZE_EQ, " G4tgrVolume::AddVisibility");

  //---------- set G4tgrPlace 
  theVisibility = G4tgrUtils::GetBool( wl[2] );
}


//-------------------------------------------------------------
void G4tgrVolume::AddRGBColour( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_GE, " G4tgrVolume::AddRGBColour");
  
  //  theRGBColour = new double[3];
  theRGBColour[0] = G4tgrUtils::GetDouble( wl[2] );
  theRGBColour[1] = G4tgrUtils::GetDouble( wl[3] );
  theRGBColour[2] = G4tgrUtils::GetDouble( wl[4] );
  if( wl.size() == 6 )
  {
    theRGBColour[3] = G4tgrUtils::GetDouble( wl[5] );
  }
}
