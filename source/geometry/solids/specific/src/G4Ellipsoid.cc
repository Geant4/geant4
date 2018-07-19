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
// $Id: G4Ellipsoid.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// class G4Ellipsoid
//
// Implementation for G4Ellipsoid class
//
// History:
//
// 10.11.99 G.Horton-Smith: first writing, based on G4Sphere class
// 25.02.05 G.Guerrieri: Modified for future Geant4 release
// 26.10.16 E.Tcherniaev: reimplemented CalculateExtent() using
//                        G4BoundingEnvelope, removed CreateRotatedVertices()
//
// --------------------------------------------------------------------

#include "globals.hh"

#include "G4Ellipsoid.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "meshdefs.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

///////////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

G4Ellipsoid::G4Ellipsoid(const G4String& pName,
                               G4double pxSemiAxis,
                               G4double pySemiAxis,
                               G4double pzSemiAxis,
                               G4double pzBottomCut,
                               G4double pzTopCut)
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0),
    fCubicVolume(0.), fSurfaceArea(0.), zBottomCut(0.), zTopCut(0.)
{
  // note: for users that want to use the full ellipsoid it is useful
  // to include a default for the cuts 

  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();

  halfCarTolerance = kCarTolerance*0.5;
  halfRadTolerance = kRadTolerance*0.5;

  // Check Semi-Axis
  if ( (pxSemiAxis<=0.) || (pySemiAxis<=0.) || (pzSemiAxis<=0.) )
  {
     std::ostringstream message;
     message << "Invalid semi-axis - " << GetName();
     G4Exception("G4Ellipsoid::G4Ellipsoid()", "GeomSolids0002",
                 FatalErrorInArgument, message);
  }
  SetSemiAxis(pxSemiAxis, pySemiAxis, pzSemiAxis);

  if ( pzBottomCut == 0 && pzTopCut == 0 )
  {
     SetZCuts(-pzSemiAxis, pzSemiAxis);
  }
  else if ( (pzBottomCut < pzSemiAxis) && (pzTopCut > -pzSemiAxis)
         && (pzBottomCut < pzTopCut) )
  {
     SetZCuts(pzBottomCut, pzTopCut);
  }
  else
  {
     std::ostringstream message;
     message << "Invalid z-coordinate for cutting plane - " << GetName();
     G4Exception("G4Ellipsoid::G4Ellipsoid()", "GeomSolids0002",
                 FatalErrorInArgument, message);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Ellipsoid::G4Ellipsoid( __void__& a )
  : G4VSolid(a), fRebuildPolyhedron(false), fpPolyhedron(0), kRadTolerance(0.),
    halfCarTolerance(0.), halfRadTolerance(0.), fCubicVolume(0.),
    fSurfaceArea(0.), xSemiAxis(0.), ySemiAxis(0.), zSemiAxis(0.),
    semiAxisMax(0.), zBottomCut(0.), zTopCut(0.)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Ellipsoid::~G4Ellipsoid()
{
  delete fpPolyhedron; fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Ellipsoid::G4Ellipsoid(const G4Ellipsoid& rhs)
  : G4VSolid(rhs),
    fRebuildPolyhedron(false), fpPolyhedron(0),
    kRadTolerance(rhs.kRadTolerance),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    xSemiAxis(rhs.xSemiAxis), ySemiAxis(rhs.ySemiAxis),
    zSemiAxis(rhs.zSemiAxis), semiAxisMax(rhs.semiAxisMax),
    zBottomCut(rhs.zBottomCut), zTopCut(rhs.zTopCut)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Ellipsoid& G4Ellipsoid::operator = (const G4Ellipsoid& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   kRadTolerance = rhs.kRadTolerance;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   fCubicVolume = rhs.fCubicVolume; fSurfaceArea = rhs.fSurfaceArea;
   xSemiAxis = rhs.xSemiAxis; ySemiAxis = rhs.ySemiAxis;
   zSemiAxis = rhs.zSemiAxis; semiAxisMax = rhs.semiAxisMax;
   zBottomCut = rhs.zBottomCut; zTopCut = rhs.zTopCut;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = 0;

   return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Ellipsoid::ComputeDimensions(G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

///////////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Ellipsoid::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double dx = GetSemiAxisMax(0);
  G4double dy = GetSemiAxisMax(1);
  G4double dz = GetSemiAxisMax(2);
  G4double zmin = std::max(-dz,GetZBottomCut());
  G4double zmax = std::min( dz,GetZTopCut());
  pMin.set(-dx,-dy,zmin);
  pMax.set( dx, dy,zmax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Ellipsoid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4Ellipsoid::CalculateExtent(const EAxis pAxis,
                             const G4VoxelLimits& pVoxelLimit,
                             const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

///////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies `in', or returns as approprate

EInside G4Ellipsoid::Inside(const G4ThreeVector& p) const
{
  G4double rad2oo,  // outside surface outer tolerance
           rad2oi;  // outside surface inner tolerance
  EInside in;

  // check this side of z cut first, because that's fast
  //
  if (p.z() < zBottomCut-halfRadTolerance) { return in=kOutside; }
  if (p.z() > zTopCut+halfRadTolerance)    { return in=kOutside; }

  rad2oo= sqr(p.x()/(xSemiAxis+halfRadTolerance))
        + sqr(p.y()/(ySemiAxis+halfRadTolerance))
        + sqr(p.z()/(zSemiAxis+halfRadTolerance));

  if (rad2oo > 1.0)  { return in=kOutside; }
    
  rad2oi= sqr(p.x()*(1.0+halfRadTolerance/xSemiAxis)/xSemiAxis)
      + sqr(p.y()*(1.0+halfRadTolerance/ySemiAxis)/ySemiAxis)
      + sqr(p.z()*(1.0+halfRadTolerance/zSemiAxis)/zSemiAxis);

  // Check radial surfaces
  //  sets `in' (already checked for rad2oo > 1.0)
  //
  if (rad2oi < 1.0)
  {
    in = ( (p.z() < zBottomCut+halfRadTolerance)
        || (p.z() > zTopCut-halfRadTolerance) ) ? kSurface : kInside;
    if ( rad2oi > 1.0-halfRadTolerance )  { in=kSurface; }
  }
  else 
  {
    in = kSurface;
  }
  return in;

}

///////////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p not protected against p=0

G4ThreeVector G4Ellipsoid::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double distR, distZBottom, distZTop;

  // normal vector with special magnitude:  parallel to normal, units 1/length
  // norm*p == 1.0 if on surface, >1.0 if outside, <1.0 if inside
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));
  G4double radius = 1.0/norm.mag();

  // approximate distance to curved surface
  //
  distR = std::fabs( (p*norm - 1.0) * radius ) / 2.0;

  // Distance to z-cut plane
  //
  distZBottom = std::fabs( p.z() - zBottomCut );
  distZTop = std::fabs( p.z() - zTopCut );

  if ( (distZBottom < distR) || (distZTop < distR) )
  {
    return G4ThreeVector(0.,0.,(distZBottom < distZTop) ? -1.0 : 1.0);
  }
  return ( norm *= radius );
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//

G4double G4Ellipsoid::DistanceToIn( const G4ThreeVector& p,
                                    const G4ThreeVector& v  ) const
{
  G4double distMin = std::min(xSemiAxis,ySemiAxis);
  const G4double dRmax = 100.*std::min(distMin,zSemiAxis);
  distMin= kInfinity;

  // check to see if Z plane is relevant
  if (p.z() <= zBottomCut+halfCarTolerance)
  {
    if (v.z() <= 0.0) { return distMin; }
    G4double distZ = (zBottomCut - p.z()) / v.z();

    if ( (distZ > -halfRadTolerance) && (Inside(p+distZ*v) != kOutside) )
    {
      // early exit since can't intercept curved surface if we reach here
      if ( std::fabs(distZ) < halfRadTolerance ) { distZ=0.; }
      return distMin= distZ;
    }
  }
  if (p.z() >= zTopCut-halfCarTolerance)
  {
    if (v.z() >= 0.0) { return distMin;}
    G4double distZ = (zTopCut - p.z()) / v.z();
    if ( (distZ > -halfRadTolerance) && (Inside(p+distZ*v) != kOutside) )
    {
      // early exit since can't intercept curved surface if we reach here
      if ( std::fabs(distZ) < halfRadTolerance ) { distZ=0.; }
      return distMin= distZ;
    }
  }
  // if fZCut1 <= p.z() <= fZCut2, then must hit curved surface

  // now check curved surface intercept
  G4double A,B,C;

  A= sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) + sqr(v.z()/zSemiAxis);
  C= sqr(p.x()/xSemiAxis) + sqr(p.y()/ySemiAxis) + sqr(p.z()/zSemiAxis) - 1.0;
  B= 2.0 * ( p.x()*v.x()/(xSemiAxis*xSemiAxis)
           + p.y()*v.y()/(ySemiAxis*ySemiAxis)
           + p.z()*v.z()/(zSemiAxis*zSemiAxis) );

  C= B*B - 4.0*A*C;
  if (C > 0.0)
  {    
    G4double distR= (-B - std::sqrt(C)) / (2.0*A);
    G4double intZ = p.z()+distR*v.z();
    if ( (distR > halfRadTolerance)
      && (intZ >= zBottomCut-halfRadTolerance)
      && (intZ <= zTopCut+halfRadTolerance) )
    { 
      distMin = distR;
    }
    else if( (distR >- halfRadTolerance)
            && (intZ >= zBottomCut-halfRadTolerance)
            && (intZ <= zTopCut+halfRadTolerance) )
    {
      // p is on the curved surface, DistanceToIn returns 0 or kInfinity:
      // DistanceToIn returns 0, if second root is positive (means going inside)
      // If second root is negative, DistanceToIn returns kInfinity (outside)
      //
      distR = (-B + std::sqrt(C) ) / (2.0*A);
      if(distR>0.) { distMin=0.; }
    }
    else
    {
      distR= (-B + std::sqrt(C)) / (2.0*A);
      intZ = p.z()+distR*v.z();
      if ( (distR > halfRadTolerance)
        && (intZ >= zBottomCut-halfRadTolerance)
        && (intZ <= zTopCut+halfRadTolerance) )
      {
        G4ThreeVector norm=SurfaceNormal(p);
        if (norm.dot(v)<0.) { distMin = distR; }
      }
    }
    if ( (distMin!=kInfinity) && (distMin>dRmax) ) 
    {                    // Avoid rounding errors due to precision issues on
                         // 64 bits systems. Split long distances and recompute
      G4double fTerm = distMin-std::fmod(distMin,dRmax);
      distMin = fTerm + DistanceToIn(p+fTerm*v,v);
    }
  }
  
  if (std::fabs(distMin)<halfRadTolerance) { distMin=0.; }
  return distMin;
} 

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Return 0 if point inside

G4double G4Ellipsoid::DistanceToIn(const G4ThreeVector& p) const
{
  G4double distR, distZ;

  // normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));
  G4double radius= 1.0/norm.mag();

  // approximate distance to curved surface ( <= actual distance )
  //
  distR= (p*norm - 1.0) * radius / 2.0;

  // Distance to z-cut plane
  //
  distZ= zBottomCut - p.z();
  if (distZ < 0.0)
  {
    distZ = p.z() - zTopCut;
  }

  // Distance to closest surface from outside
  //
  if (distZ < 0.0)
  {
    return (distR < 0.0) ? 0.0 : distR;
  }
  else if (distR < 0.0)
  {
    return distZ;
  }
  else
  {
    return (distZ < distR) ? distZ : distR;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p,
                                    const G4ThreeVector& v,
                                    const G4bool calcNorm,
                                          G4bool *validNorm,
                                          G4ThreeVector *n  ) const
{
  G4double distMin;
  enum surface_e {kPlaneSurf, kCurvedSurf, kNoSurf} surface;
  
  distMin= kInfinity;
  surface= kNoSurf;

  // check to see if Z plane is relevant
  //
  if (v.z() < 0.0)
  {
    G4double distZ = (zBottomCut - p.z()) / v.z();
    if (distZ < 0.0)
    {
      distZ= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    distMin= distZ;
    surface= kPlaneSurf;
  }
  if (v.z() > 0.0)
  {
    G4double distZ = (zTopCut - p.z()) / v.z();
    if (distZ < 0.0)
    {
      distZ= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    distMin= distZ;
    surface= kPlaneSurf;
  }

  // normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector nearnorm(p.x()/(xSemiAxis*xSemiAxis),
                         p.y()/(ySemiAxis*ySemiAxis),
                         p.z()/(zSemiAxis*zSemiAxis));
  
  // now check curved surface intercept
  //
  G4double A,B,C;
  
  A= sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) + sqr(v.z()/zSemiAxis);
  C= (p * nearnorm) - 1.0;
  B= 2.0 * (v * nearnorm);

  C= B*B - 4.0*A*C;
  if (C > 0.0)
  {
    G4double distR= (-B + std::sqrt(C) ) / (2.0*A);
    if (distR < 0.0)
    {
      distR= 0.0;
      if (!calcNorm) {return 0.0;}
    }
    if (distR < distMin)
    {
      distMin= distR;
      surface= kCurvedSurf;
    }
  }

  // set normal if requested
  //
  if (calcNorm)
  {
    if (surface == kNoSurf)
    {
      *validNorm = false;
    }
    else
    {
      *validNorm = true;
      switch (surface)
      {
        case kPlaneSurf:
          *n= G4ThreeVector(0.,0.,(v.z() > 0.0 ? 1. : -1.));
          break;
        case kCurvedSurf:
        {
          G4ThreeVector pexit= p + distMin*v;
          G4ThreeVector truenorm(pexit.x()/(xSemiAxis*xSemiAxis),
                                 pexit.y()/(ySemiAxis*ySemiAxis),
                                 pexit.z()/(zSemiAxis*zSemiAxis));
          truenorm *= 1.0/truenorm.mag();
          *n= truenorm;
        } break;
        default:           // Should never reach this case ...
          DumpInfo();
          std::ostringstream message;
          G4int oldprc = message.precision(16);
          message << "Undefined side for valid surface normal to solid."
                  << G4endl
                  << "Position:"  << G4endl
                  << "   p.x() = "   << p.x()/mm << " mm" << G4endl
                  << "   p.y() = "   << p.y()/mm << " mm" << G4endl
                  << "   p.z() = "   << p.z()/mm << " mm" << G4endl
                  << "Direction:" << G4endl << G4endl
                  << "   v.x() = "   << v.x() << G4endl
                  << "   v.y() = "   << v.y() << G4endl
                  << "   v.z() = "   << v.z() << G4endl
                  << "Proposed distance :" << G4endl
                  << "   distMin = "    << distMin/mm << " mm";
          message.precision(oldprc);
          G4Exception("G4Ellipsoid::DistanceToOut(p,v,..)",
                      "GeomSolids1002", JustWarning, message);
          break;
      }
    }
  }
   
  return distMin;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Ellipsoid::DistanceToOut(const G4ThreeVector& p) const
{
  G4double distR, distZ;

#ifdef G4SPECSDEBUG
  if( Inside(p) == kOutside )
  {
     DumpInfo();
     std::ostringstream message;
     G4int oldprc = message.precision(16);
     message << "Point p is outside !?" << G4endl
             << "Position:"  << G4endl
             << "   p.x() = "   << p.x()/mm << " mm" << G4endl
             << "   p.y() = "   << p.y()/mm << " mm" << G4endl
             << "   p.z() = "   << p.z()/mm << " mm";
     message.precision(oldprc) ;
     G4Exception("G4Ellipsoid::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, message);
  }
#endif

  // Normal vector:  parallel to normal, magnitude 1/(characteristic radius)
  //
  G4ThreeVector norm(p.x()/(xSemiAxis*xSemiAxis),
                     p.y()/(ySemiAxis*ySemiAxis),
                     p.z()/(zSemiAxis*zSemiAxis));

  // the following is a safe inlined "radius= min(1.0/norm.mag(),p.mag())
  //
  G4double radius= p.mag();
  G4double tmp= norm.mag();
  if ( (tmp > 0.0) && (1.0 < radius*tmp) ) {radius = 1.0/tmp;}

  // Approximate distance to curved surface ( <= actual distance )
  //
  distR = (1.0 - p*norm) * radius / 2.0;
    
  // Distance to z-cut plane
  //
  distZ = p.z() - zBottomCut;
  if (distZ < 0.0) {distZ= zTopCut - p.z();}

  // Distance to closest surface from inside
  //
  if ( (distZ < 0.0) || (distR < 0.0) )
  {
    return 0.0;
  }
  else
  {
    return (distZ < distR) ? distZ : distR;
  }
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType

G4GeometryType G4Ellipsoid::GetEntityType() const
{
  return G4String("G4Ellipsoid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Ellipsoid::Clone() const
{
  return new G4Ellipsoid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Ellipsoid::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Ellipsoid\n"
     << " Parameters: \n"

     << "    semi-axis x: " << xSemiAxis/mm << " mm \n"
     << "    semi-axis y: " << ySemiAxis/mm << " mm \n"
     << "    semi-axis z: " << zSemiAxis/mm << " mm \n"
     << "    max semi-axis: " << semiAxisMax/mm << " mm \n"
     << "    lower cut plane level z: " << zBottomCut/mm << " mm \n"
     << "    upper cut plane level z: " << zTopCut/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Ellipsoid::GetPointOnSurface() const
{
  G4double aTop, aBottom, aCurved, chose, xRand, yRand, zRand, phi;
  G4double cosphi, sinphi, costheta, sintheta, alpha, beta, max1, max2, max3;

  max1  = xSemiAxis > ySemiAxis ? xSemiAxis : ySemiAxis;
  max1  = max1 > zSemiAxis ? max1 : zSemiAxis;
  if (max1 == xSemiAxis)      { max2 = ySemiAxis; max3 = zSemiAxis; }
  else if (max1 == ySemiAxis) { max2 = xSemiAxis; max3 = zSemiAxis; }
  else                        { max2 = xSemiAxis; max3 = ySemiAxis; }

  phi   = G4RandFlat::shoot(0.,twopi);
  
  cosphi = std::cos(phi);   sinphi = std::sin(phi);
  costheta = G4RandFlat::shoot(zBottomCut,zTopCut)/zSemiAxis;
  sintheta = std::sqrt(1.-sqr(costheta));
  
  alpha = 1.-sqr(max2/max1); beta  = 1.-sqr(max3/max1);
  
  aTop    = pi*xSemiAxis*ySemiAxis*(1 - sqr(zTopCut/zSemiAxis));
  aBottom = pi*xSemiAxis*ySemiAxis*(1 - sqr(zBottomCut/zSemiAxis));
  
  // approximation
  // from:" http://www.citr.auckland.ac.nz/techreports/2004/CITR-TR-139.pdf"
  aCurved = 4.*pi*max1*max2*(1.-1./6.*(alpha+beta)-
                            1./120.*(3.*sqr(alpha)+2.*alpha*beta+3.*sqr(beta)));

  aCurved *= 0.5*(1.2*zTopCut/zSemiAxis - 1.2*zBottomCut/zSemiAxis);
  
  if( ( zTopCut >= zSemiAxis && zBottomCut <= -1.*zSemiAxis )
   || ( zTopCut == 0 && zBottomCut ==0 ) )
  {
    aTop = 0; aBottom = 0;
  }
  
  chose = G4RandFlat::shoot(0.,aTop + aBottom + aCurved); 
  
  if(chose < aCurved)
  { 
    xRand = xSemiAxis*sintheta*cosphi;
    yRand = ySemiAxis*sintheta*sinphi;
    zRand = zSemiAxis*costheta;
    return G4ThreeVector (xRand,yRand,zRand); 
  }
  else if(chose >= aCurved && chose < aCurved + aTop)
  {
    xRand = G4RandFlat::shoot(-1.,1.)*xSemiAxis
          * std::sqrt(1-sqr(zTopCut/zSemiAxis));
    yRand = G4RandFlat::shoot(-1.,1.)*ySemiAxis
          * std::sqrt(1.-sqr(zTopCut/zSemiAxis)-sqr(xRand/xSemiAxis));
    zRand = zTopCut;
    return G4ThreeVector (xRand,yRand,zRand);
  }
  else
  {
    xRand = G4RandFlat::shoot(-1.,1.)*xSemiAxis
          * std::sqrt(1-sqr(zBottomCut/zSemiAxis));
    yRand = G4RandFlat::shoot(-1.,1.)*ySemiAxis
          * std::sqrt(1.-sqr(zBottomCut/zSemiAxis)-sqr(xRand/xSemiAxis)); 
    zRand = zBottomCut;
    return G4ThreeVector (xRand,yRand,zRand);
  }
}

/////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Ellipsoid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4VisExtent G4Ellipsoid::GetExtent() const
{
  // Define the sides of the box into which the G4Ellipsoid instance would fit.
  //
  return G4VisExtent (-semiAxisMax, semiAxisMax,
                      -semiAxisMax, semiAxisMax,
                      -semiAxisMax, semiAxisMax);
}

G4Polyhedron* G4Ellipsoid::CreatePolyhedron () const
{
  return new G4PolyhedronEllipsoid(xSemiAxis, ySemiAxis, zSemiAxis,
                                   zBottomCut, zTopCut);
}

G4Polyhedron* G4Ellipsoid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      G4AutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }
  return fpPolyhedron;
}
