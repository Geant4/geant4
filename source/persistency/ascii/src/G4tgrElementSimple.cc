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
// $Id: G4tgrElementSimple.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrElementSimple

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrElementSimple.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgrElementSimple::G4tgrElementSimple()
  : theZ(0.), theA(0.)
{
}


// -------------------------------------------------------------------------
G4tgrElementSimple::~G4tgrElementSimple()
{
}


// -------------------------------------------------------------------------
G4tgrElementSimple::G4tgrElementSimple( const std::vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ,
                           "G4tgrElementSimple::G4tgrElementSimple");

  theType = "ElementSimple";
  theName = G4tgrUtils::GetString( wl[1] );
  theSymbol = G4tgrUtils::GetString( wl[2] );
  theZ = G4tgrUtils::GetInt( wl[3] );
  theA = G4tgrUtils::GetDouble( wl[4], g/mole);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrElementSimple& obj)
{
  os << "G4tgrElementSimple= " << obj.theName
     << " Z = " << obj.theZ << " A= " << obj.theA << G4endl;

  return os;
}
