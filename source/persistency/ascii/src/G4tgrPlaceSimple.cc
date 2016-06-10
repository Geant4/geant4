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
// $Id: G4tgrPlaceSimple.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrPlaceSimple

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrPlaceSimple.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrVolume.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgrPlaceSimple::G4tgrPlaceSimple()
{
}


// -------------------------------------------------------------------------
G4tgrPlaceSimple::~G4tgrPlaceSimple()
{
}


// -------------------------------------------------------------------------
G4tgrPlaceSimple::G4tgrPlaceSimple( const std::vector<G4String>& wl )
{
  theType = "PlaceSimple";

  G4int wl7 = -1;
  if( wl.size() == 8 )  // for assembly volume placement,
  {                     // there is no copy number
    //---------- set the copy number
    theCopyNo = G4tgrUtils::GetInt( wl[2] );
    wl7 = 0;
  }

  //---------- set the parent name
  theParentName = G4tgrUtils::GetString( wl[3+wl7] ); 

  //---------- set the position with respect to parent
  thePlace = G4ThreeVector( G4tgrUtils::GetDouble(wl[5+wl7])*mm,
                            G4tgrUtils::GetDouble(wl[6+wl7])*mm,
                            G4tgrUtils::GetDouble(wl[7+wl7])*mm ); 

  //---------- set the rotation matrix name
  theRotMatName = G4tgrUtils::GetString(wl[4+wl7]);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrPlaceSimple& obj)
{
  os << "G4tgrPlaceSimple=  in " << obj.theParentName
     << " Position= " << obj.thePlace
     << " RotMatName= " << obj.theRotMatName << G4endl;
     
  return os;
}
