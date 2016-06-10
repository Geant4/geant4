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
// $Id: G4tgrVolumeAssembly.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrVolumeAssembly

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrVolumeAssembly.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolumeMgr.hh"
#include "G4tgrPlaceSimple.hh"
#include "G4tgrPlaceDivRep.hh"
#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrFileReader.hh"
#include "G4tgrMessenger.hh"

//-------------------------------------------------------------
G4tgrVolumeAssembly::G4tgrVolumeAssembly()
{
}


//-------------------------------------------------------------
G4tgrVolumeAssembly::~G4tgrVolumeAssembly()
{
}


//-------------------------------------------------------------
G4tgrVolumeAssembly::G4tgrVolumeAssembly( const std::vector<G4String>& wl ) 
{
  theType = "VOLAssembly";

  //---------- set name 
  theName = G4tgrUtils::GetString( wl[1] ); 

  G4int nVol = G4tgrUtils::GetInt( wl[2] ); 

  G4tgrUtils::CheckWLsize( wl, 3+nVol*5, WLSIZE_GE,
                           "G4tgrVolumeAssembly::G4tgrVolumeAssembly" );


  for(G4int ii=0; ii<nVol*5; ii+=5)
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 2 )
    {
         G4cout << " G4tgrVolumeAssembly::G4tgrVolumeAssembly() -"
                << " Adding component: " << ii << " - " << wl[ii+3] << G4endl;
    }
#endif
    theComponentNames.push_back(G4tgrUtils::GetString( wl[3+ii+0]));
    theComponentRMs.push_back(G4tgrUtils::GetString( wl[3+ii+1]));
    theComponentPos.push_back(G4ThreeVector(G4tgrUtils::GetDouble(wl[3+ii+2]),
                                            G4tgrUtils::GetDouble(wl[3+ii+3]),
                                            G4tgrUtils::GetDouble(wl[3+ii+4])));
  }
  theVisibility = 1;
  theRGBColour = new G4double[4];
  for (size_t ii=0; ii<4; ii++)  { theRGBColour[ii] = -1.; }

  theSolid = 0;

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif

}


//-------------------------------------------------------------------------
G4tgrPlace* G4tgrVolumeAssembly::AddPlace( const std::vector<G4String>& wl )
{
  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_EQ, " G4tgrVolumeAssembly::AddPlace");

  //---------- set G4tgrPlace 
  G4tgrPlaceSimple* pl = new G4tgrPlaceSimple( wl );

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


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrVolumeAssembly& obj)
{
  os << "G4tgrVolumeAssembly= " << obj.theName;

  for( size_t ii = 0; ii < obj.theComponentNames.size(); ii++ )
  {
    os << obj.theComponentNames[ii] << " RotMatName= "
       << obj.theComponentRMs[ii] << " Position= "
       << obj.theComponentPos[ii].x() << " "
       << obj.theComponentPos[ii].y() << " "
       << obj.theComponentPos[ii].z();
  }
  os << G4endl;

  return os;
}
