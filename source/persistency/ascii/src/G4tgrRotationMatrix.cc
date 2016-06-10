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
// $Id: G4tgrRotationMatrix.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrRotationMatrix

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrRotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrRotationMatrixFactory.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgrRotationMatrix::G4tgrRotationMatrix()
  : theName("Rotation-Matrix"), theInputType(rm9)
{
}


// -------------------------------------------------------------------------
G4tgrRotationMatrix::~G4tgrRotationMatrix()
{
}


// -------------------------------------------------------------------------
G4tgrRotationMatrix::G4tgrRotationMatrix( const std::vector<G4String>& wl )
  : theInputType(rm9) 
{
  theName = G4tgrUtils::GetString( wl[1] );

  switch( wl.size() )
  {
    case 5:
      theInputType = rm3;
      break;
    case 8:
      theInputType = rm6;
      break;
    case 11:
      theInputType = rm9;
      break;
    default:
      G4Exception("G4tgrRotationMatrix::G4tgrRotationMatrix()",
                  "InvalidMatrix", FatalException,
                  "Input line must have 5, 8 or 11 words.");
      break;
  }
 
  //-------- Fill matrix values
  size_t siz = wl.size() - 2;
  for( size_t ii = 0; ii < siz; ii++)
  {
    if( siz == 9 )
    {
      theValues.push_back( G4tgrUtils::GetDouble( wl[ii+2] ) );
    }
    else
    {
      theValues.push_back( G4tgrUtils::GetDouble( wl[ii+2] , deg ) );
    }
  }
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << " G4tgrRotationMatrix::G4tgrRotationMatrix() - Created: "
           << theName << G4endl;
    for( size_t ii = 0; ii < siz; ii++)
    {
      G4cout << " " << theValues[ii];
    }
    G4cout << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrRotationMatrix& obj)
{
  os << "G4tgrRotationMatrix= " << obj.theName
     << " InputTyep = " << obj.theInputType << " VALUES= ";

  for( size_t ii = 0; ii < obj.theValues.size(); ii++ )
  {
    os << obj.theValues[ii] << " ";
  }

  os << G4endl;

  return os;
}
