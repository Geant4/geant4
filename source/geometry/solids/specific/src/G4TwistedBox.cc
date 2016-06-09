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
// $Id: G4TwistedBox.cc,v 1.11 2005/12/06 09:22:13 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedBox.cc
//
// Author:
//
//   10/11/2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4TwistedBox.hh"
#include "G4Polyhedron.hh"

G4TwistedBox::G4TwistedBox( const G4String& pName,
                                  G4double  pPhiTwist,
                                  G4double  pDx, 
                                  G4double  pDy, 
                                  G4double  pDz )
  : G4VTwistedFaceted( pName, pPhiTwist,pDz,0.,0.,
                       pDy, pDx, pDx, pDy, pDx, pDx,0. )
{
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4TwistedBox::G4TwistedBox( __void__& a )
  : G4VTwistedFaceted(a)
{
}

G4TwistedBox::~G4TwistedBox()
{
}

std::ostream& G4TwistedBox::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedBox\n"
     << " Parameters: \n"
     << "    pDx = "  << GetXHalfLength()/cm << " cm" << G4endl
     << "    pDy = "  << GetYHalfLength()/cm << " cm" << G4endl
     << "    pDz = "  << GetZHalfLength()/cm << " cm" << G4endl
     << "    pPhiTwist = " << GetPhiTwist()/degree << " deg" << G4endl 
     << "-----------------------------------------------------------\n";

  return os;
}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedBox::GetEntityType() const
{
  return G4String("G4TwistedBox");
}
