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
// $Id: G4TwistedTrd.cc,v 1.2 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
//   Oliver Link (Oliver.Link@cern.ch)
//
// History:
//
//   18/03/2005 - Created.
//
// --------------------------------------------------------------------

#include "G4TwistedTrd.hh"
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

G4TwistedTrd::~G4TwistedTrd()
{
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
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4TwistedTrd::CreatePolyhedron () const 
{
  // Returns simple Trd for now !!
  //
  return new G4PolyhedronTrd2 (GetDx1(),GetDx3(),GetDy1(),GetDy2(),GetDz());
}

//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTrd::GetEntityType() const
{
  return G4String("G4TwistedTrd");
}
