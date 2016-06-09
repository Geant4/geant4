//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TwistedTrap.cc,v 1.13 2005/12/06 09:22:13 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
