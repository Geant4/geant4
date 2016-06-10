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
// $Id: G4tgrPlaceParameterisation.cc 73052 2013-08-15 15:45:36Z gcosmo $
//
//
// class G4tgrPlaceParameterisation

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrPlaceParameterisation.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrMessenger.hh"

//-------------------------------------------------------------
G4tgrPlaceParameterisation::G4tgrPlaceParameterisation()
{
}


//-------------------------------------------------------------
G4tgrPlaceParameterisation::~G4tgrPlaceParameterisation()
{
}


//-------------------------------------------------------------
G4tgrPlaceParameterisation::
G4tgrPlaceParameterisation( const std::vector<G4String>& wl )
{
  theType = "PlaceParam";

  //---------- Check for exact number of words read 
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_GE,
                           "G4tgrPlaceParameterisation::ConstructVolume" );
  
  //---------- the copy No
  theCopyNo = G4tgrUtils::GetInt( wl[2] )-1;

  //---------- set the parent name
  theParentName = G4tgrUtils::GetString( wl[3] ); 

  //---------- set the type
  theParamType = G4tgrUtils::GetString( wl[4] );

  //---------- set the rotation matrix name
  theRotMatName = G4tgrUtils::GetString(wl[5]);
 
  //---------- set the extra data 
  for( size_t ii = 6; ii < wl.size(); ii++)
  {
    theExtraData.push_back( G4tgrUtils::GetDouble(wl[ii]) );
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os,
                         const G4tgrPlaceParameterisation& obj)
{
  os << "G4tgrPlaceParameterisation= in " << obj.theParentName
     << " ParamTyep = " << obj.theParamType
     << " RotMatName= " << obj.theRotMatName << " EXTRA DATA= ";
  for( size_t ii = 0; ii < obj.theExtraData.size(); ii++ )
  {
    os << obj.theExtraData[ii] << " ";
  }
  os << G4endl;

  return os;
}

