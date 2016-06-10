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
// $Id: G4tgrMaterialMixture.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgrMaterialMixture

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgrMaterialMixture.hh"

#include "G4SystemOfUnits.hh"
#include "G4tgrUtils.hh"
#include "G4tgrMessenger.hh"

//-----------------------------------------------------------
G4tgrMaterialMixture::G4tgrMaterialMixture()
{
}


//-----------------------------------------------------------
G4tgrMaterialMixture::~G4tgrMaterialMixture()
{
}


//-----------------------------------------------------------
G4tgrMaterialMixture::G4tgrMaterialMixture(const G4String& matType,
                                           const std::vector<G4String>& wl)
{
  //---------- Check for miminum number of words read 
  G4tgrUtils::CheckWLsize( wl, 6, WLSIZE_GE,
                           "G4tgrMaterialMixture::G4tgrMaterialMixture" );

  theMateType = matType;
  
  //---------- Fill private data 
  theName = G4tgrUtils::GetString( wl[1] );
  theDensity = std::fabs(G4tgrUtils::GetDouble( wl[2], g/cm3 ) );
  theNoComponents = G4tgrUtils::GetInt( wl[3] );

  G4tgrUtils::CheckWLsize( wl, 4+theNoComponents*2, WLSIZE_GE,
                           "G4tgrMaterialMixture::G4tgrMaterialMixture" );
  for(G4int ii=0; ii<theNoComponents; ii++)
  {
#ifdef G4VERBOSE
    if( G4tgrMessenger::GetVerboseLevel() >= 3 )
    {
         G4cout << " G4tgrMaterialMixture::G4tgrMaterialMixture() -"
                << " adding component: " << wl[ii*2+4] << " Fraction= "
                << G4tgrUtils::GetDouble(wl[ii*2+1+4]) << G4endl;
    }
#endif
    theComponents.push_back(  G4tgrUtils::GetString( wl[ii*2+4] ) );
    theFractions.push_back( G4tgrUtils::GetDouble(wl[ii*2+1+4]) );
  }

#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 1 )
  {
     G4cout << " Created " << *this << G4endl;
  }
#endif
}


//-----------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgrMaterialMixture& mate) 
{
  os << "G4tgrMaterialMixture=: " << mate.theName << G4endl
     << "density= " << mate.theDensity/g*cm3
     << " g/cm3. Number of Components: " << mate.theNoComponents << G4endl;
  for (G4int ii=0; ii<mate.theNoComponents; ii++)
  {
    os << '\t' << mate.theComponents[ii]
       << '\t' << mate.theFractions[ii] << G4endl;
  }
  return os;
}
