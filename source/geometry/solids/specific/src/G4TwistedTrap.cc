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
// $Id: G4TwistedTrap.cc 79492 2014-03-05 15:25:30Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrap.cc
//
// Author:
//
//   10/11/2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4TwistedTrap.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedron.hh"

G4TwistedTrap::G4TwistedTrap( const G4String &pName,
                                    G4double  pPhiTwist,
                                    G4double  pDx1,
                                    G4double  pDx2,
                                    G4double  pDy,
                                    G4double  pDz)
  : G4VTwistedFaceted( pName, pPhiTwist,pDz,0.,0.,
                       pDy, pDx1, pDx2, pDy, pDx1, pDx2,0. )
{
}

G4TwistedTrap::
G4TwistedTrap(const G4String &pName,      // Name of instance
                    G4double  pPhiTwist,  // twist angle
                    G4double  pDz,        // half z length
                    G4double  pTheta,  // direction between end planes
                    G4double  pPhi,    // defined by polar and azimuthal angles
                    G4double  pDy1,    // half y length at -pDz
                    G4double  pDx1,    // half x length at -pDz,-pDy
                    G4double  pDx2,    // half x length at -pDz,+pDy
                    G4double  pDy2,    // half y length at +pDz
                    G4double  pDx3,    // half x length at +pDz,-pDy
                    G4double  pDx4,    // half x length at +pDz,+pDy
                    G4double  pAlph    // tilt angle
              )
  : G4VTwistedFaceted( pName, pPhiTwist, pDz, pTheta,
                       pPhi, pDy1, pDx1, pDx2, pDy2, pDx3, pDx4, pAlph )
{
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4TwistedTrap::G4TwistedTrap( __void__& a )
  : G4VTwistedFaceted(a)
{
}

G4TwistedTrap::~G4TwistedTrap()
{
}

// Copy constructor
//
G4TwistedTrap::G4TwistedTrap(const G4TwistedTrap& rhs)
  : G4VTwistedFaceted(rhs)
{
  fpPolyhedron = GetPolyhedron();
}

// Assignment operator
//
G4TwistedTrap& G4TwistedTrap::operator = (const G4TwistedTrap& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VTwistedFaceted::operator=(rhs);
   fpPolyhedron = GetPolyhedron();

   return *this;
}

std::ostream& G4TwistedTrap::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedTrap\n"
     << " Parameters: \n"
     << "    Twist angle         = " << GetPhiTwist()/degree << " deg"
     << G4endl 
     << "    Polar Angle Theta   = " << GetPolarAngleTheta()/degree << " deg"
     << G4endl 
     << "    Azimuthal Angle Phi = " << GetAzimuthalAnglePhi()/degree << " deg"
     << G4endl 
     << "    pDy1 = " << GetY1HalfLength()/cm << " cm" << G4endl
     << "    pDx1 = " << GetX1HalfLength()/cm << " cm" << G4endl
     << "    pDx2 = " << GetX2HalfLength()/cm << " cm" << G4endl
     << "    pDy2 = " << GetY2HalfLength()/cm << " cm" << G4endl
     << "    pDx3 = " << GetX3HalfLength()/cm << " cm" << G4endl
     << "    pDx4 = " << GetX4HalfLength()/cm << " cm" << G4endl
     << "    pDz = "  << GetZHalfLength()/cm << " cm" << G4endl
     << "    Tilt Angle Alpha    = " << GetTiltAngleAlpha()/degree << " deg"
     << G4endl 
     << "-----------------------------------------------------------\n";

  return os;
}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTrap::GetEntityType() const
{
  return G4String("G4TwistedTrap");
}

//=====================================================================
//* Clone -------------------------------------------------------------

G4VSolid* G4TwistedTrap::Clone() const
{
  return new G4TwistedTrap(*this);
}
