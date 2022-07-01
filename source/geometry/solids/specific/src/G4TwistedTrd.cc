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
// G4TwistedTrd
//
// Author: 18/03/2005 - O.Link (Oliver.Link@cern.ch)
// --------------------------------------------------------------------

#include "G4TwistedTrd.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedron.hh"

//=====================================================================
//* Constructor -------------------------------------------------------

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

//=====================================================================
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4TwistedTrd::G4TwistedTrd( __void__& a )
  : G4VTwistedFaceted(a)
{
}

//=====================================================================
//* Destructor --------------------------------------------------------

G4TwistedTrd::~G4TwistedTrd()
{
}

//=====================================================================
//* Copy constructor --------------------------------------------------

G4TwistedTrd::G4TwistedTrd(const G4TwistedTrd& rhs)
  : G4VTwistedFaceted(rhs)
{
  fpPolyhedron = GetPolyhedron();
}

//=====================================================================
//* Assignment operator -----------------------------------------------

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

//=====================================================================
//* StreamInfo --------------------------------------------------------

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

//=====================================================================
//* GetCubicVolume ----------------------------------------------------

double G4TwistedTrd::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    G4double x1 = GetX1HalfLength();
    G4double x2 = GetX2HalfLength();
    G4double y1 = GetY1HalfLength();
    G4double y2 = GetY2HalfLength();
    G4double h = 2.*GetZHalfLength();
    fCubicVolume = h*((x1 + x2)*(y1 + y2) + (x2 - x1)*(y2 - y1)/3.);
  }
  return fCubicVolume;
}

//=====================================================================
//* GetSurfaceArea ----------------------------------------------------

double G4TwistedTrd::GetSurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    G4double ang = GetPhiTwist();
    G4double x1 = GetX1HalfLength();
    G4double x2 = GetX2HalfLength();
    G4double y1 = GetY1HalfLength();
    G4double y2 = GetY2HalfLength();
    G4double h = 2.*GetZHalfLength();
    G4double hh = h*h;
    G4double delX = x2 - x1;
    G4double delY = y2 - y1;
    if (ang == 0.)
    {
      G4double hx = std::sqrt(delY*delY + hh);
      G4double hy = std::sqrt(delX*delX + hh);
      return fSurfaceArea =
        2.*(x1 + x2)*hx + 2.*(y1 + y2)*hy + 4.*(x1*y1 + x2*y2);
    }

    // compute area of x-faces
    G4double U1, U2, V1, V2;
    G4double areaX = 0.;
    U1 = delY + x1*ang;
    U2 = delY + x2*ang;
    V1 = delY - x1*ang;
    V2 = delY - x2*ang;
    if (std::abs(delX) < kCarTolerance) // case x1 == x2
    {
      areaX = (U1*std::sqrt(hh + U1*U1) + hh*std::asinh(U1/h) -
               V1*std::sqrt(hh + V1*V1) - hh*std::asinh(V1/h))/ang;
    }
    else
    {
      // U contribution
      areaX += ((hh + U2*U2)*std::sqrt(hh + U2*U2) -
                (hh + U1*U1)*std::sqrt(hh + U1*U1))/3.
        + hh*(U2*std::asinh(U2/h) - U1*std::asinh(U1/h))
        - hh*(std::sqrt(hh + U2*U2) - std::sqrt(hh + U1*U1));
      // V contribution
      areaX += ((hh + V2*V2)*std::sqrt(hh + V2*V2) -
                (hh + V1*V1)*std::sqrt(hh + V1*V1))/3.
        + hh*(V2*std::asinh(V2/h) - V1*std::asinh(V1/h))
        - hh*(std::sqrt(hh + V2*V2) - std::sqrt(hh + V1*V1));
      areaX /= delX*ang*ang;
    }

    // compute area of y-faces
    G4double areaY = 0.;
    U1 = delX + y1*ang;
    U2 = delX + y2*ang;
    V1 = delX - y1*ang;
    V2 = delX - y2*ang;
    if (std::abs(delY) < kCarTolerance) // case y1 == y2
    {
      areaY = (U1*std::sqrt(hh + U1*U1) + hh*std::asinh(U1/h) -
               V1*std::sqrt(hh + V1*V1) - hh*std::asinh(V1/h))/ang;
    }
    else
    {
      // U contribution
      areaY += ((hh + U2*U2)*std::sqrt(hh + U2*U2) -
                (hh + U1*U1)*std::sqrt(hh + U1*U1))/3.
        + hh*(U2*std::asinh(U2/h) - U1*std::asinh(U1/h))
        - hh*(std::sqrt(hh + U2*U2) - std::sqrt(hh + U1*U1));
      // V contribution
      areaY += ((hh + V2*V2)*std::sqrt(hh + V2*V2) -
                (hh + V1*V1)*std::sqrt(hh + V1*V1))/3.
        + hh*(V2*std::asinh(V2/h) - V1*std::asinh(V1/h))
        - hh*(std::sqrt(hh + V2*V2) - std::sqrt(hh + V1*V1));
      areaY /= delY*ang*ang;
    }
    fSurfaceArea = areaX + areaY + 4.*(x1*y1 + x2*y2);
  }
  return fSurfaceArea;
}
