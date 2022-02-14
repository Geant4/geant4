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
  G4ThreeVector pt[8];
  CheckParameters();
  GetVertices(pt);
  CheckPlanarity(pt);
}

G4UTrap::G4UTrap( const G4String& pName,
                  const G4ThreeVector pt[8] )
  : G4USolid(pName, new UTrap(pName))
{
  // Start with check of centering - the center of gravity trap line
  // should cross the origin of frame
  if (!(   pt[0].z() < 0
        && pt[0].z() == pt[1].z()
        && pt[0].z() == pt[2].z()
        && pt[0].z() == pt[3].z()

        && pt[4].z() > 0
        && pt[4].z() == pt[5].z()
        && pt[4].z() == pt[6].z()
        && pt[4].z() == pt[7].z()

        && std::abs( pt[0].z() + pt[4].z() ) < kCarTolerance

        && pt[0].y() == pt[1].y()
        && pt[2].y() == pt[3].y()
        && pt[4].y() == pt[5].y()
        && pt[6].y() == pt[7].y()

        && std::abs(pt[0].y()+pt[2].y()+pt[4].y()+pt[6].y()) < kCarTolerance
        && std::abs(pt[0].x()+pt[1].x()+pt[4].x()+pt[5].x() +
                    pt[2].x()+pt[3].x()+pt[6].x()+pt[7].x()) < kCarTolerance ))
  {
    std::ostringstream message;
    message << "Invalid vertice coordinates for Solid: " << GetName();
    G4Exception("G4UTrap::G4UTrap()", "GeomSolids0002",
                FatalException, message);
  }

  SetPlanes(pt);
  CheckParameters();
  CheckPlanarity(pt);
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pZ,
                        G4double pY,
                        G4double pX, G4double pLTX )
  : G4USolid(pName, new UTrap(pName, pZ, pY, pX, pLTX))
{
  CheckParameters();
}

G4UTrap::G4UTrap( const G4String& pName,
                        G4double pdx1,  G4double pdx2,
                        G4double pdy1,  G4double pdy2,
                        G4double pdz )
  : G4USolid(pName, new UTrap(pName, pdx1, pdx2, pdy1, pdy2, pdz))
{
  CheckParameters();
}

G4UTrap::G4UTrap(const G4String& pName,
                       G4double pdx, G4double pdy, G4double pdz,
                       G4double pAlpha, G4double pTheta, G4double pPhi )
  : G4USolid(pName, new UTrap(pName, pdx, pdy, pdz, pAlpha, pTheta, pPhi))
{
  CheckParameters();
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
// Accessors
//
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
G4double G4UTrap::GetTanAlpha1() const
{
  return Base_t::GetTanAlpha1();
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
  return GetDx4();
}
G4double G4UTrap::GetTanAlpha2() const
{
  return Base_t::GetTanAlpha2();
}
G4double G4UTrap::GetPhi() const       
{
  return Base_t::GetPhi();
}
G4double G4UTrap::GetTheta() const
{
  return Base_t::GetTheta();
}
G4double G4UTrap::GetAlpha1() const
{
  return Base_t::GetAlpha1();
}
G4double G4UTrap::GetAlpha2() const
{
  return Base_t::GetAlpha2();
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

//////////////////////////////////////////////////////////////////////////
//
// Modifier
//
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

  G4ThreeVector pt[8];
  CheckParameters();
  GetVertices(pt);
  CheckPlanarity(pt);
}

/////////////////////////////////////////////////////////////////////////
//
// Set parameters using eight vertices
//
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
// Check dimensions
//
void G4UTrap::CheckParameters() const
{
  G4double fDz  = GetZHalfLength();
  G4double fDy1 = GetYHalfLength1();
  G4double fDx1 = GetXHalfLength1();
  G4double fDx2 = GetXHalfLength2();
  G4double fDy2 = GetYHalfLength2();
  G4double fDx3 = GetXHalfLength3();
  G4double fDx4 = GetXHalfLength4();

  if (fDz<=0 ||
      fDy1<=0 || fDx1<=0 || fDx2<=0 ||
      fDy2<=0 || fDx3<=0 || fDx4<=0)
  {
    std::ostringstream message;
    message << "Invalid Length Parameters for Solid: " << GetName()
            << "\n  X - " <<fDx1<<", "<<fDx2<<", "<<fDx3<<", "<<fDx4
            << "\n  Y - " <<fDy1<<", "<<fDy2
            << "\n  Z - " <<fDz;
    G4Exception("G4UTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Compute coordinates of vertices
//
void G4UTrap::GetVertices(G4ThreeVector pt[8]) const
{
  G4double fDz      = GetZHalfLength();
  G4double fDy1     = GetYHalfLength1();
  G4double fDx1     = GetXHalfLength1();
  G4double fDx2     = GetXHalfLength2();
  G4double fDy2     = GetYHalfLength2();
  G4double fDx3     = GetXHalfLength3();
  G4double fDx4     = GetXHalfLength4();
  G4double phi      = GetPhi();
  G4double theta    = GetTheta();
  G4double fTalpha1 = GetTanAlpha1();
  G4double fTalpha2 = GetTanAlpha2();

  G4double DzTthetaCphi = fDz*std::tan(theta)*std::cos(phi);
  G4double DzTthetaSphi = fDz*std::tan(theta)*std::sin(phi);
  G4double Dy1Talpha1   = fDy1*fTalpha1;
  G4double Dy2Talpha2   = fDy2*fTalpha2;

  pt[0].set(-DzTthetaCphi-Dy1Talpha1-fDx1,-DzTthetaSphi-fDy1,-fDz);
  pt[1].set(-DzTthetaCphi-Dy1Talpha1+fDx1,-DzTthetaSphi-fDy1,-fDz);
  pt[2].set(-DzTthetaCphi+Dy1Talpha1-fDx2,-DzTthetaSphi+fDy1,-fDz);
  pt[3].set(-DzTthetaCphi+Dy1Talpha1+fDx2,-DzTthetaSphi+fDy1,-fDz);
  pt[4].set( DzTthetaCphi-Dy2Talpha2-fDx3, DzTthetaSphi-fDy2, fDz);
  pt[5].set( DzTthetaCphi-Dy2Talpha2+fDx3, DzTthetaSphi-fDy2, fDz);
  pt[6].set( DzTthetaCphi+Dy2Talpha2-fDx4, DzTthetaSphi+fDy2, fDz);
  pt[7].set( DzTthetaCphi+Dy2Talpha2+fDx4, DzTthetaSphi+fDy2, fDz);
}

/////////////////////////////////////////////////////////////////////////
//
// Check planarity of lateral planes
//
void G4UTrap::CheckPlanarity(const G4ThreeVector pt[8]) const
{
  constexpr G4int iface[4][4] = { {0,4,5,1}, {2,3,7,6}, {0,2,6,4}, {1,5,7,3} };
  const static G4String side[4] = { "~-Y", "~+Y", "~-X", "~+X" };

  for (G4int i=0; i<4; ++i)
  {
    TrapSidePlane plane = GetSidePlane(i);
    G4double dmax = 0;
    for (G4int k=0; k<4; ++k)
    {
      const G4ThreeVector p = pt[iface[i][k]];
      G4double dist = plane.a*p.x() + plane.b*p.y() + plane.c*p.z() + plane.d;
      if (std::abs(dist) > std::abs(dmax)) dmax = dist;
    }
    if (std::abs(dmax) > 1000 * kCarTolerance)
    {
      std::ostringstream message;
      message << "Side face " << side[i] << " is not planar for solid: "
              << GetName() << "\nDiscrepancy: " << dmax/mm << " mm\n";
      StreamInfo(message);
      G4Exception("G4UTrap::CheckPlanarity()", "GeomSolids0002",
                  FatalException, message);
    }
  }
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
  return new G4PolyhedronTrap(GetZHalfLength(), GetTheta(), GetPhi(),
                              GetYHalfLength1(),
                              GetXHalfLength1(), GetXHalfLength2(), GetAlpha1(),
                              GetYHalfLength2(),
                              GetXHalfLength3(), GetXHalfLength4(), GetAlpha2());
}

#endif  // G4GEOM_USE_USOLIDS
