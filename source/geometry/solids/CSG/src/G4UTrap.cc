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
// Implementation for G4UTrap wrapper class
//
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#include "G4Trap.hh"
#include "G4UTrap.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4VPVParameterisation.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdz,
                        G4double pTheta, G4double pPhi,
                        G4double pdy1, G4double pdx1, G4double pdx2,
                        G4double pAlp1,
                        G4double pdy2, G4double pdx3, G4double pdx4,
                        G4double pAlp2 )
  : G4USolid(pName, new UTrap(pName, pdz, pTheta, pPhi,
                              pdy1, pdx1, pdx2, pAlp1, pdy2, pdx3, pdx4, pAlp2))
{
}

G4UTrap::G4UTrap( const G4String& pName,
                  const G4ThreeVector pt[8] )
  : G4USolid(pName, new UTrap(pName))
{
  SetPlanes(pt);
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pZ,
                        G4double pY,
                        G4double pX, G4double pLTX )
  : G4USolid(pName, new UTrap(pName, pZ, pY, pX, pLTX))
{
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdx1,  G4double pdx2,
                        G4double pdy1,  G4double pdy2,
                        G4double pdz )
  : G4USolid(pName, new UTrap(pName, pdx1, pdx2, pdy1, pdy2, pdz))
{
}

G4UTrap::G4UTrap(const G4String& pName,
                       G4double pdx, G4double pdy, G4double pdz,
                       G4double pAlpha, G4double pTheta, G4double pPhi )
  : G4USolid(pName, new UTrap(pName, pdx, pdy, pdz, pAlpha, pTheta, pPhi))
{
}

G4UTrap::G4UTrap( const G4String& pName )
  : G4USolid(pName, new UTrap(pName))
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTrap::G4UTrap( __void__& a )
  : G4USolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTrap::~G4UTrap()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTrap::G4UTrap(const G4UTrap& rhs)
  : G4USolid(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTrap& G4UTrap::operator = (const G4UTrap& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4USolid::operator=(rhs);

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers

G4double G4UTrap::GetZHalfLength() const
{
  return GetShape()->GetZHalfLength();
}
G4double G4UTrap::GetYHalfLength1() const
{
  return GetShape()->GetYHalfLength1();
}
G4double G4UTrap::GetXHalfLength1() const
{
  return GetShape()->GetXHalfLength1();
}
G4double G4UTrap::GetXHalfLength2() const
{
  return GetShape()->GetXHalfLength2();
}
G4double G4UTrap::GetTanAlpha1() const
{
  return GetShape()->GetTanAlpha1();
}
G4double G4UTrap::GetYHalfLength2() const
{
  return GetShape()->GetYHalfLength2();
}
G4double G4UTrap::GetXHalfLength3() const
{
  return GetShape()->GetXHalfLength3();
}
G4double G4UTrap::GetXHalfLength4() const
{
  return GetShape()->GetXHalfLength4();
}
G4double G4UTrap::GetTanAlpha2() const
{
  return GetShape()->GetTanAlpha2();
}
TrapSidePlane G4UTrap::GetSidePlane(G4int n) const
{
  UTrapSidePlane iplane = GetShape()->GetSidePlane(n);
  TrapSidePlane oplane = {iplane.a, iplane.b, iplane.c, iplane.d };
  return oplane;
}
G4ThreeVector G4UTrap::GetSymAxis() const
{
  UVector3 axis = GetShape()->GetSymAxis();
  return G4ThreeVector(axis.x(), axis.y(), axis.z());
}

void G4UTrap::SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                               G4double pDy1, G4double pDx1, G4double pDx2,
                               G4double pAlp1,
                               G4double pDy2, G4double pDx3, G4double pDx4,
                               G4double pAlp2)
{
  GetShape()->SetAllParameters(pDz, pTheta, pPhi,
                               pDy1, pDx1, pDx2, pAlp1,
                               pDy2, pDx3, pDx4, pAlp2);
  fRebuildPolyhedron = true;
}

void G4UTrap::SetPlanes(const G4ThreeVector pt[8])
{
  UVector3 upt[8];
  for (unsigned int i=0; i<8; ++i)
  {
    upt[i] = UVector3(pt[i].x(), pt[i].y(), pt[i].z());
  }
  GetShape()->SetPlanes(upt);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UTrap::ComputeDimensions(      G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Trap*)this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4UTrap::Clone() const
{
  return new G4UTrap(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization
//
G4Polyhedron* G4UTrap::CreatePolyhedron() const
{
  G4double fTthetaSphi = GetShape()->GetThetaSphi();
  G4double fTthetaCphi = GetShape()->GetThetaCphi();
  G4double phi = std::atan2(fTthetaSphi, fTthetaCphi);
  G4double alpha1 = std::atan(GetTanAlpha1());
  G4double alpha2 = std::atan(GetTanAlpha2());
  G4double theta = std::atan(std::sqrt(fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi));

  return new G4PolyhedronTrap(GetZHalfLength(), theta, phi,
                              GetYHalfLength1(),
                              GetXHalfLength1(), GetXHalfLength2(), alpha1,
                              GetYHalfLength2(),
                              GetXHalfLength3(), GetXHalfLength4(), alpha2);
}

#endif  // G4GEOM_USE_USOLIDS
