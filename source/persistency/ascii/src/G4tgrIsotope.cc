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
// $Id: G4tgrIsotope.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrIsotope

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrIsotope.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"


//-------------------------------------------------------------
G4tgrIsotope::G4tgrIsotope()
  : theName(""), theZ(0), theN(0), theA(0.)
{
}


//-------------------------------------------------------------
G4tgrIsotope::~G4tgrIsotope()
{
}


//-------------------------------------------------------------
G4tgrIsotope::G4tgrIsotope( const std::vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ, "G4tgrIsotope::G4tgIstotope");

  theName = G4tgrUtils::GetString( wl[1] );
  theZ = G4tgrUtils::GetInt( wl[2] );
  theN = G4tgrUtils::GetInt( wl[3] );
  theA = G4tgrUtils::GetDouble( wl[4], g/mole);

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif
}

// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrIsotope& obj)
{
  os << "G4tgrIsotope= " << obj.theName
     << " Z = " << obj.theZ
     << " N= " << obj.theN
     << " A= " << obj.theA << G4endl;

  return os;
}
