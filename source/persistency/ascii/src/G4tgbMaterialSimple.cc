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
// $Id: G4tgbMaterialSimple.cc 66363 2012-12-18 09:12:54Z gcosmo $
//
//
// class G4tgbMaterialSimple

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#include "G4tgbMaterialSimple.hh"

#include "G4PhysicalConstants.hh"
#include "G4tgrMaterialSimple.hh"
#include "G4tgrMessenger.hh"

// -------------------------------------------------------------------------
G4tgbMaterialSimple::G4tgbMaterialSimple()
  : theZ(0.), theA(0.)
{
}


// -------------------------------------------------------------------------
G4tgbMaterialSimple::~G4tgbMaterialSimple()
{
}


// -------------------------------------------------------------------------
G4tgbMaterialSimple::G4tgbMaterialSimple( G4tgrMaterial* hgmate)
{
  theTgrMate = hgmate;
  G4tgrMaterialSimple* matesimp = static_cast<G4tgrMaterialSimple*>(hgmate);
  theZ = matesimp->GetZ();
  theA = matesimp->GetA();
}


// -------------------------------------------------------------------------
G4Material* G4tgbMaterialSimple::BuildG4Material()
{
  //----- construct new G4Material with no components (only itself)

  G4Material* mate = new G4Material( GetName(), GetZ(), GetA(),
                                     theTgrMate->GetDensity(),
                                     kStateUndefined, STP_Temperature );
  
#ifdef G4VERBOSE
  if( G4tgrMessenger::GetVerboseLevel() >= 2 )
  {
    G4cout << "  Constructing new G4Material simple: " << *mate << G4endl;
  }
#endif

  return mate;
}


// -------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const G4tgbMaterialSimple& mate) 
{
  os << "Simple Material: " << mate.GetName() << G4endl
     << " Z = " << mate.GetZ() 
     << " A = " << mate.GetA() 
     << " density = " << mate.GetDensity() << G4endl;
  return os;
}
