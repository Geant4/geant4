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
// $Id: G4tgrPlaceDivRep.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrPlaceDivRep

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrPlaceDivRep.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrVolume.hh"

//-------------------------------------------------------------
G4tgrPlaceDivRep::G4tgrPlaceDivRep()
  : theNDiv(0), theWidth(0.), theAxis(kUndefined),
    theOffset(0.), theDivType(DivByNdivAndWidth) 
{
}


//-------------------------------------------------------------
G4tgrPlaceDivRep::~G4tgrPlaceDivRep()
{
}

//-------------------------------------------------------------
G4tgrPlaceDivRep::G4tgrPlaceDivRep( const std::vector<G4String>& wl ) 
{
  theDivType = DivByNdivAndWidth;

  // Name parent axis nrep width offset
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE,
                           "G4tgrPlaceDivRep::G4tgrPlaceDivRep" );
  G4tgrUtils::CheckWLsize( wl, 7, WLSIZE_LE,
                           "G4tgrPlaceDivRep::G4tgrPlaceDivRep" );

  theParentName = G4tgrUtils::GetString(wl[2]); 
  theAxis = BuildAxis( G4tgrUtils::GetString(wl[3]) ); 
  theNDiv = G4tgrUtils::GetInt( wl[4] );
  theWidth = G4tgrUtils::GetDouble(wl[5])*mm; // check if it is deg
  if( wl.size() == 7 )
  {
    theOffset = G4tgrUtils::GetDouble(wl[6])*mm;  // check if it is deg
  }
  else
  {
    theOffset = 0.;
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif

}


//-------------------------------------------------------------
EAxis G4tgrPlaceDivRep::BuildAxis( const G4String& axisName ) 
{
  if( axisName == "X" ) {
    return kXAxis;
  } else if( axisName == "Y" ) {
    return kYAxis;
  } else if( axisName == "Z" ) {
    return kZAxis;
  } else if( axisName == "R" ) {
    return kRho;
  } else if( axisName == "PHI" ) {
    return kPhi;
  }
  else
  {
    G4String ErrMessage = "Axis type not found: " + axisName
                        + ". Only valid axis are: X, Y, Z, R, PHI !";
    G4Exception("G4tgrVolumeDivision::GetReplicaAxis()",
                "InvalidAxis", FatalException, ErrMessage);
  }
  return kXAxis; // to avoid warning errors  
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrPlaceDivRep& obj)
{
  os << "G4tgrPlaceDivRep= in " << obj.theParentName
     << " NDiv= " << obj.theNDiv << " Width= " << obj.theWidth
     << " Axis= " << obj.theAxis << " Offset= " << obj.theOffset
     << " DivType= " << obj.theDivType << G4endl;

  return os;
}
