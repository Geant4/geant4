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
// $Id: G4tgrMaterialSimple.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrMaterialSimple

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrMaterialSimple.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"
#include "G4UIcommand.hh"

//-------------------------------------------------------------
G4tgrMaterialSimple::G4tgrMaterialSimple()
  : name("MaterialSimple"), theA(0.), theZ(0.)
{
}


//-------------------------------------------------------------
G4tgrMaterialSimple::~G4tgrMaterialSimple()
{
}


//-------------------------------------------------------------
G4tgrMaterialSimple::G4tgrMaterialSimple(const G4String& matType,
                                         const std::vector<G4String>& wl)
  : name("MaterialSimple")
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 5, WLSIZE_EQ,
                           "G4tgrMaterialSimple::G4tgrMaterialSimple");

  theMateType = matType;

  //---------- Fill private data 
  theName = G4tgrUtils::GetString( wl[1] );
  theZ = G4tgrUtils::GetDouble( wl[2] );
  theA = G4tgrUtils::GetDouble( wl[3], g/mole);
  theDensity = G4tgrUtils::GetDouble( wl[4], g/cm3);
  theNoComponents = 0;
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif
}


//-------------------------------------------------------------
const G4String& G4tgrMaterialSimple::GetComponent(G4int ii) const
{ 
  G4String ErrMessage = "Should never be called for a MaterialSimple - i:"
                      + G4UIcommand::ConvertToString(ii);
  G4Exception("G4tgrMaterialSimple::GetComponent()",
              "InvalidCall", FatalException, ErrMessage);

  return name;  // dummy, to avoid compilation warnings...
}


//-------------------------------------------------------------
G4double G4tgrMaterialSimple::GetFraction(G4int ii)
{
  G4String ErrMessage = "Should never be called for a MaterialSimple - i:"
                      + G4UIcommand::ConvertToString(ii);
  G4Exception("G4tgrMaterialSimple::GetFraction()",
              "InvalidCall", FatalException, ErrMessage);

  return 0.; // dummy, to avoid compilation warnings...
}


//-------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrMaterialSimple& mate) 
{
  os << "G4tgrMaterialSimple= " << mate.theName 
     << " Z " << mate.theZ << " A " << mate.theA 
     << "density= " << mate.theDensity/g*cm3
     << " g/cm3. Number of Components: " << mate.theNoComponents << G4endl;
  return os;
}
