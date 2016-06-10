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
// $Id: G4TwistedTrd.cc 79492 2014-03-05 15:25:30Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrd.cc
//
// Author:
//
//   18/03/2005 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4TwistedTrd.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedron.hh"

G4TwistedTrd::G4TwistedTrd( const G4String& pName,
                                  G4double  pDx1, 
                                  G4double  pDx2,
                                  G4double  pDy1, 
                                  G4double  pDy2,
                                  G4double  pDz,
                                  G4double  pPhiTwist )
  : G4VTwistedFaceted( pName, pPhiTwist,pDz,0.,0.,
                       pDy1, pDx1, pDx1, pDy2, pDx2, pDx2,0.)
{
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4TwistedTrd::G4TwistedTrd( __void__& a )
  : G4VTwistedFaceted(a)
{
}

G4TwistedTrd::~G4TwistedTrd()
{
}

// Copy constructor
//
G4TwistedTrd::G4TwistedTrd(const G4TwistedTrd& rhs)
  : G4VTwistedFaceted(rhs)
{
  fpPolyhedron = GetPolyhedron();
}

// Assignment operator
//
G4TwistedTrd& G4TwistedTrd::operator = (const G4TwistedTrd& rhs) 
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

std::ostream& G4TwistedTrd::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedTrd\n"
     << " Parameters: \n"
     << "    pDx1 = " << GetX1HalfLength()/cm << " cm" << G4endl
     << "    pDx2 = " << GetX2HalfLength()/cm << " cm" << G4endl
     << "    pDy1 = " << GetY1HalfLength()/cm << " cm" << G4endl
     << "    pDy2 = " << GetY2HalfLength()/cm << " cm" << G4endl
     << "    pDz = "  << GetZHalfLength()/cm << " cm" << G4endl
     << "    pPhiTwist = " << GetPhiTwist()/degree << " deg" << G4endl 
     << "-----------------------------------------------------------\n";

  return os;
}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTrd::GetEntityType() const
{
  return G4String("G4TwistedTrd");
}

//=====================================================================
//* Clone -------------------------------------------------------------

G4VSolid* G4TwistedTrd::Clone() const
{
  return new G4TwistedTrd(*this);
}
