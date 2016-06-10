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
// $Id: G4tgrVolume.cc 66363 2012-12-18 09:12:54Z gcosmo $
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
#include "G4UIcommand.hh"

//-------------------------------------------------------------
G4tgrVolume::G4tgrVolume()
  : theName(""), theType(""),
    theMaterialName(""), theSolid(0), theVisibility(false),
    theRGBColour(0), theCheckOverlaps(false)
{
}


//-------------------------------------------------------------
G4tgrVolume::~G4tgrVolume()
{
  delete [] theRGBColour;
}


//-------------------------------------------------------------
G4tgrVolume::G4tgrVolume( const std::vector<G4String>& wl) 
{
  theType = "VOLSimple";

  //---------- set name 
  theName = G4tgrUtils::GetString( wl[1] ); 

  theVisibility = 1;
  theRGBColour = new G4double[4];
  for(size_t ii=0; ii<4; ii++)  { theRGBColour[ii] = -1.; }
  theCheckOverlaps = 0;

  if( wl.size() != 4 )
  {
    //:VOLU tag to build a volume creating solid and material
    //---------- set material name
    theMaterialName = G4tgrUtils::GetString( wl[wl.size()-1] );
    
    //---------- create only vector<double> of theSolidParams
    theSolid = G4tgrVolumeMgr::GetInstance()->CreateSolid( wl, 1 );

#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 1 )
      {
        G4cout << "Created from new solid: " 
               << *this << G4endl;
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
      if( G4tgrMessenger::GetVerboseLevel() >= 1 )
        {
          G4cout << "Created from existing solid: " 
                 << *this << G4endl;
        }
#endif
    }
}


//-------------------------------------------------------------------------
G4tgrVolume::G4tgrVolume( const G4tgrVolume& vol )
{
  theName = vol.GetName();   
  theType = vol.GetType();
  theMaterialName = vol.GetMaterialName();   
  theSolid = vol.GetSolid();
  thePlacements  = vol.GetPlacements();
  theVisibility   = vol.GetVisibility();
  theRGBColour   = vol.GetRGBColour();
  theCheckOverlaps = vol.GetCheckOverlaps();
}


//-------------------------------------------------------------------------
G4tgrVolume* G4tgrVolume::GetVolume( G4int ii ) const
{
  G4String ErrMessage = "Should only be called for composite solids... "
                      + G4UIcommand::ConvertToString(ii);
  G4Exception("G4tgrVolume::GetVolume()", "InvalidCall",
              FatalException, ErrMessage);
  return 0;
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
    G4cout << " G4tgrVolume:  New placement: " << thePlacements.size()
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

  if( (wl.size() == 7) && (G4tgrUtils::GetDouble(wl[6]) != 0.)
                       && (wl[3] != "PHI") )
  {
    G4Exception("G4tgrVolume::AddPlaceReplica",
                "Offset set for replica not along PHI, it will not be used",
                JustWarning,
                G4String("Volume "+wl[1]+" in volume "+wl[2]).c_str());
  }
  
  //---------- set G4tgrPlace 
  G4tgrPlaceDivRep* pl = new G4tgrPlaceDivRep( wl );
  pl->SetType("PlaceReplica");
  pl->SetVolume( this );
  thePlacements.push_back( pl ); 

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrVolume:  New placement replica: " << thePlacements.size()
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
    G4cout << " G4tgrVolume:  New placement Param: " << thePlacements.size()
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

  //---------- Set visibility
  theVisibility = G4tgrUtils::GetBool( wl[2] );
}


//-------------------------------------------------------------
void G4tgrVolume::AddRGBColour( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_GE, " G4tgrVolume::AddRGBColour");
  
  //---------- Set RGB colour
  theRGBColour[0] = G4tgrUtils::GetDouble( wl[2] );
  theRGBColour[1] = G4tgrUtils::GetDouble( wl[3] );
  theRGBColour[2] = G4tgrUtils::GetDouble( wl[4] );
  ///--------- Set transparency
  if( wl.size() == 6 )
  {
    theRGBColour[3] = G4tgrUtils::GetDouble( wl[5] );
  }
}


//-------------------------------------------------------------
void G4tgrVolume::AddCheckOverlaps( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 3, WLSIZE_GE, " G4tgrVolume::AddCheckOverlaps");
  
  ///--------- Set check overlaps 
  theCheckOverlaps = G4tgrUtils::GetBool( wl[2] );

}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrVolume& obj)
{
  os << "G4tgrVolume= " << obj.theName << " Type= " << obj.theType
     << " Material= " << obj.theMaterialName
     << " Visibility " << obj.theVisibility 
     << " Colour " << (obj.theRGBColour)[0] << " "
                   << (obj.theRGBColour)[1] << " "
                   << (obj.theRGBColour)[2] << " "
                   << (obj.theRGBColour)[3] << " "
     << " CheckOverlaps " << obj.theCheckOverlaps
     << " N placements " << obj.thePlacements.size() << G4endl;
      
  return os;
}
