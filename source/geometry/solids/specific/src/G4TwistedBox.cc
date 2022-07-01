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
// G4TwistedBox implementation
//
// Author:  10/11/2004 - O.Link (Oliver.Link@cern.ch)
// --------------------------------------------------------------------

#include "G4TwistedBox.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedron.hh"

//=====================================================================
//* Constructor -------------------------------------------------------

G4TwistedBox::G4TwistedBox( const G4String& pName,
                                  G4double  pPhiTwist,
                                  G4double  pDx,
                                  G4double  pDy,
                                  G4double  pDz )
  : G4VTwistedFaceted( pName, pPhiTwist,pDz,0.,0.,
                       pDy, pDx, pDx, pDy, pDx, pDx,0. )
{
}

//=====================================================================
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4TwistedBox::G4TwistedBox( __void__& a )
  : G4VTwistedFaceted(a)
{
}

//=====================================================================
//* Destructor --------------------------------------------------------

G4TwistedBox::~G4TwistedBox()
{
}

//=====================================================================
//* Copy constructor --------------------------------------------------

G4TwistedBox::G4TwistedBox(const G4TwistedBox& rhs)
  : G4VTwistedFaceted(rhs)
{
  fpPolyhedron = GetPolyhedron();
}

//=====================================================================
//* Assignment operator -----------------------------------------------

G4TwistedBox& G4TwistedBox::operator = (const G4TwistedBox& rhs)
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

//=====================================================================
//* StreamInfo --------------------------------------------------------

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

//=====================================================================
//* Clone -------------------------------------------------------------

G4VSolid* G4TwistedBox::Clone() const
{
  return new G4TwistedBox(*this);
}

//=====================================================================
//* GetCubicVolume ----------------------------------------------------

double G4TwistedBox::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    fCubicVolume = 8.*GetXHalfLength()*GetYHalfLength()*GetZHalfLength();
  }
  return fCubicVolume;
}

//=====================================================================
//* GetSurfaceArea ----------------------------------------------------

double G4TwistedBox::GetSurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    G4double ang = GetPhiTwist();
    G4double dx = GetXHalfLength();
    G4double dy = GetYHalfLength();
    G4double dz = GetZHalfLength();
    if (ang == 0.)
    {
      fSurfaceArea = 8.*(dx*dy + dx*dz + dy*dz);
    }
    else
    {
      G4double h = 2.*dz;
      G4double hh = h*h;
      G4double dxang = dx*ang;
      G4double dyang = dy*ang;
      fSurfaceArea = 8.*dx*dy +
        2.*(dx*std::sqrt(hh + dxang*dxang) + hh*std::asinh(dxang/h)/ang) +
        2.*(dy*std::sqrt(hh + dyang*dyang) + hh*std::asinh(dyang/h)/ang);
    }
  }
  return fSurfaceArea;
}
