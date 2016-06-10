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
// $Id: G4tgrElementFromIsotopes.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrElementFromIsotopes

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrElementFromIsotopes.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"


// -------------------------------------------------------------------------
G4tgrElementFromIsotopes::G4tgrElementFromIsotopes()
  : theNoIsotopes(0)
{
}

// -------------------------------------------------------------------------
G4tgrElementFromIsotopes::~G4tgrElementFromIsotopes()
{
}

// -------------------------------------------------------------------------
G4tgrElementFromIsotopes::
G4tgrElementFromIsotopes( const std::vector<G4String>& wl ) 
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE,
                          "G4tgrElementFromIsotopes::G4tgrElementFromIsotopes");
    //:ELEM_FROM_ISOT NAME SYMBOL N_ISOT (ISOT_NAME ISOT_ABUNDANCE)

  theType = "ElementFromIsotopes";
  theName = G4tgrUtils::GetString( wl[1] );
  theSymbol = G4tgrUtils::GetString( wl[2] );
  theNoIsotopes = G4tgrUtils::GetInt( wl[3] );

  for( G4int ii = 0; ii < theNoIsotopes; ii++ )
  {
    theComponents.push_back( G4tgrUtils::GetString( wl[4+ii*2] ) ); 
    theAbundances.push_back( G4tgrUtils::GetDouble( wl[4+ii*2+1] ) );
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
    G4cout << " Created " << *this << G4endl;
  }
#endif
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrElementFromIsotopes& obj)
{
  os << "G4tgrElementFromIsotopes= " << obj.theName
     << " N isotopes " << obj.theNoIsotopes
     << " COMPONENTS " << G4endl;
  for(size_t ii = 0; ii < obj.theComponents.size(); ii++ )
  { 
    os << obj.theComponents[ii] << " : " << obj.theAbundances[ii];
  }
  os << G4endl;

  return os;
}
