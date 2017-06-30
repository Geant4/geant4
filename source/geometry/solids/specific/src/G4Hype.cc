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
// $Id: G4Hype.cc 104316 2017-05-24 13:04:23Z gcosmo $
// $Original: G4Hype.cc,v 1.0 1998/06/09 16:57:50 safai Exp $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4Hype.cc
//
// --------------------------------------------------------------------
//
// Authors: 
//      Ernesto Lamanna (Ernesto.Lamanna@roma1.infn.it) &
//      Francesco Safai Tehrani (Francesco.SafaiTehrani@roma1.infn.it)
//      Rome, INFN & University of Rome "La Sapienza",  9 June 1998.
//
// --------------------------------------------------------------------

#include "G4Hype.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4ClippablePolygon.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include <cmath>

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

// Constructor - check parameters, and fills protected data members
G4Hype::G4Hype(const G4String& pName,
                     G4double newInnerRadius,
                     G4double newOuterRadius,
                     G4double newInnerStereo,
                     G4double newOuterStereo,
                     G4double newHalfLenZ)
  : G4VSolid(pName), fCubicVolume(0.), fSurfaceArea(0.),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fHalfTol = 0.5*kCarTolerance;

  // Check z-len
  //
  if (newHalfLenZ<=0)
  {
    std::ostringstream message;
    message << "Invalid Z half-length - " << GetName() << G4endl
            << "        Invalid Z half-length: "
            << newHalfLenZ/mm << " mm";
    G4Exception("G4Hype::G4Hype()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
  halfLenZ=newHalfLenZ;   

  // Check radii
  //
  if (newInnerRadius<0 || newOuterRadius<0)
  {
    std::ostringstream message;
    message << "Invalid radii - " << GetName() << G4endl
            << "        Invalid radii !  Inner radius: "
            << newInnerRadius/mm << " mm" << G4endl
            << "                         Outer radius: "
            << newOuterRadius/mm << " mm";
    G4Exception("G4Hype::G4Hype()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
  if (newInnerRadius >= newOuterRadius)
  {
    std::ostringstream message;
    message << "Outer > inner radius - " << GetName() << G4endl
            << "        Invalid radii !  Inner radius: "
            << newInnerRadius/mm << " mm" << G4endl
            << "                         Outer radius: "
            << newOuterRadius/mm << " mm";
    G4Exception("G4Hype::G4Hype()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  innerRadius=newInnerRadius;
  outerRadius=newOuterRadius;

  innerRadius2=innerRadius*innerRadius;
  outerRadius2=outerRadius*outerRadius;
    
  SetInnerStereo( newInnerStereo );
  SetOuterStereo( newOuterStereo );
}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Hype::G4Hype( __void__& a  )
  : G4VSolid(a), innerRadius(0.), outerRadius(0.), halfLenZ(0.), innerStereo(0.),
    outerStereo(0.), tanInnerStereo(0.), tanOuterStereo(0.), tanInnerStereo2(0.),
    tanOuterStereo2(0.), innerRadius2(0.), outerRadius2(0.), endInnerRadius2(0.),
    endOuterRadius2(0.), endInnerRadius(0.), endOuterRadius(0.),
    fCubicVolume(0.), fSurfaceArea(0.), fHalfTol(0.),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
}


//
// Destructor
//
G4Hype::~G4Hype()
{
  delete fpPolyhedron; fpPolyhedron = 0;
}


//
// Copy constructor
//
G4Hype::G4Hype(const G4Hype& rhs)
  : G4VSolid(rhs), innerRadius(rhs.innerRadius),
    outerRadius(rhs.outerRadius), halfLenZ(rhs.halfLenZ),
    innerStereo(rhs.innerStereo), outerStereo(rhs.outerStereo),
    tanInnerStereo(rhs.tanInnerStereo), tanOuterStereo(rhs.tanOuterStereo),
    tanInnerStereo2(rhs.tanInnerStereo2), tanOuterStereo2(rhs.tanOuterStereo2),
    innerRadius2(rhs.innerRadius2), outerRadius2(rhs.outerRadius2),
    endInnerRadius2(rhs.endInnerRadius2), endOuterRadius2(rhs.endOuterRadius2),
    endInnerRadius(rhs.endInnerRadius), endOuterRadius(rhs.endOuterRadius),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fHalfTol(rhs.fHalfTol), fRebuildPolyhedron(false), fpPolyhedron(0)
{
}


//
// Assignment operator
//
G4Hype& G4Hype::operator = (const G4Hype& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   innerRadius = rhs.innerRadius; outerRadius = rhs.outerRadius;
   halfLenZ = rhs.halfLenZ;
   innerStereo = rhs.innerStereo; outerStereo = rhs.outerStereo;
   tanInnerStereo = rhs.tanInnerStereo; tanOuterStereo = rhs.tanOuterStereo;
   tanInnerStereo2 = rhs.tanInnerStereo2; tanOuterStereo2 = rhs.tanOuterStereo2;
   innerRadius2 = rhs.innerRadius2; outerRadius2 = rhs.outerRadius2;
   endInnerRadius2 = rhs.endInnerRadius2; endOuterRadius2 = rhs.endOuterRadius2;
   endInnerRadius = rhs.endInnerRadius; endOuterRadius = rhs.endOuterRadius;
   fCubicVolume = rhs.fCubicVolume; fSurfaceArea = rhs.fSurfaceArea;
   fHalfTol = rhs.fHalfTol;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = 0;

   return *this;
}


//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4Hype::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Hype::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin.set(-endOuterRadius,-endOuterRadius,-halfLenZ);
  pMax.set( endOuterRadius, endOuterRadius, halfLenZ);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Hype::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Hype::CalculateExtent(const EAxis pAxis,
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


//
// Decides whether point is inside, outside or on the surface
//
EInside G4Hype::Inside(const G4ThreeVector& p) const
{
  //
  // Check z extents: are we outside?
  //
  const G4double absZ(std::fabs(p.z()));
  if (absZ > halfLenZ + fHalfTol) return kOutside;
  
  //
  // Check outer radius
  //
  const G4double oRad2(HypeOuterRadius2(absZ));
  const G4double xR2( p.x()*p.x()+p.y()*p.y() );
  
  if (xR2 > oRad2 + kCarTolerance*endOuterRadius) return kOutside;
  
  if (xR2 > oRad2 - kCarTolerance*endOuterRadius) return kSurface;
  
  if (InnerSurfaceExists())
  {
    //
    // Check inner radius
    //
    const G4double iRad2(HypeInnerRadius2(absZ));
    
    if (xR2 < iRad2 - kCarTolerance*endInnerRadius) return kOutside;
    
    if (xR2 < iRad2 + kCarTolerance*endInnerRadius) return kSurface;
  }
  
  //
  // We are inside in radius, now check endplate surface
  //
  if (absZ > halfLenZ - fHalfTol) return kSurface;
  
  return kInside;
}



//
// return the normal unit vector to the Hyperbolical Surface at a point 
// p on (or nearly on) the surface
//
G4ThreeVector G4Hype::SurfaceNormal( const G4ThreeVector& p ) const
{
  //
  // Which of the three or four surfaces are we closest to?
  //
  const G4double absZ(std::fabs(p.z()));
  const G4double distZ(absZ - halfLenZ);
  const G4double dist2Z(distZ*distZ);
  
  const G4double xR2( p.x()*p.x()+p.y()*p.y() );
  const G4double dist2Outer( std::fabs(xR2 - HypeOuterRadius2(absZ)) );
  
  if (InnerSurfaceExists())
  {
    //
    // Has inner surface: is this closest?
    //
    const G4double dist2Inner( std::fabs(xR2 - HypeInnerRadius2(absZ)) );
    if (dist2Inner < dist2Z && dist2Inner < dist2Outer)
      return G4ThreeVector( -p.x(), -p.y(), p.z()*tanInnerStereo2 ).unit();
  }

  //
  // Do the "endcaps" win?
  //
  if (dist2Z < dist2Outer) 
    return G4ThreeVector( 0.0, 0.0, p.z() < 0 ? -1.0 : 1.0 );
    
    
  //
  // Outer surface wins
  //
  return G4ThreeVector( p.x(), p.y(), -p.z()*tanOuterStereo2 ).unit();
}


//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection,
//   or intersection distance <= tolerance
//
// Calculating the intersection of a line with the surfaces
// is fairly straight forward. The difficult problem is dealing
// with the intersections of the surfaces in a consistent manner, 
// and this accounts for the complicated logic.
//
G4double G4Hype::DistanceToIn( const G4ThreeVector& p,
                               const G4ThreeVector& v ) const
{
  //
  // Quick test. Beware! This assumes v is a unit vector!
  //
  if (std::fabs(p.x()*v.y() - p.y()*v.x()) > endOuterRadius+kCarTolerance)
    return kInfinity;
  
  //
  // Take advantage of z symmetry, and reflect throught the
  // z=0 plane so that pz is always positive
  //
  G4double pz(p.z()), vz(v.z());
  if (pz < 0)
  {
    pz = -pz;
    vz = -vz;
  }

  //
  // We must be very careful if we don't want to
  // create subtle leaks at the edges where the
  // hyperbolic surfaces connect to the endplate.
  // The only reliable way to do so is to make sure
  // that the decision as to when a track passes
  // over the edge of one surface is exactly the
  // same decision as to when a track passes into the
  // other surface. By "exact", we don't mean algebraicly
  // exact, but we mean the same machine instructions
  // should be used.
  //
  G4bool couldMissOuter(true),
         couldMissInner(true),
         cantMissInnerCylinder(false);
  
  //
  // Check endplate intersection
  //
  G4double sigz = pz-halfLenZ;
  
  if (sigz > -fHalfTol)    // equivalent to: if (pz > halfLenZ - fHalfTol)
  {
    //
    // We start in front of the endplate (within roundoff)
    // Correct direction to intersect endplate?
    //
    if (vz >= 0)
    {
      //
      // Nope. As long as we are far enough away, we
      // can't intersect anything
      //
      if (sigz > 0) return kInfinity;
      
      //
      // Otherwise, we may still hit a hyperbolic surface
      // if the point is on the hyperbolic surface (within tolerance)
      //
      G4double pr2 = p.x()*p.x() + p.y()*p.y();
      if (pr2 > endOuterRadius2 + kCarTolerance*endOuterRadius)
        return kInfinity;
      
      if (InnerSurfaceExists())
      {
        if (pr2 < endInnerRadius2 - kCarTolerance*endInnerRadius)
          return kInfinity;
        if ( (pr2 < endOuterRadius2 - kCarTolerance*endOuterRadius)
          && (pr2 > endInnerRadius2 + kCarTolerance*endInnerRadius) )
          return kInfinity;
      }
      else
      {
        if (pr2 < endOuterRadius2 - kCarTolerance*endOuterRadius)
          return kInfinity;
      }
    }
    else
    {
      //
      // Where do we intersect at z = halfLenZ?
      //
      G4double q = -sigz/vz;
      G4double xi = p.x() + q*v.x(),
               yi = p.y() + q*v.y();
         
      //
      // Is this on the endplate? If so, return s, unless
      // we are on the tolerant surface, in which case return 0
      //
      G4double pr2 = xi*xi + yi*yi;
      if (pr2 <= endOuterRadius2)
      {
        if (InnerSurfaceExists())
        {
          if (pr2 >= endInnerRadius2) return (sigz < fHalfTol) ? 0 : q;
          //
          // This test is sufficient to ensure that the
          // trajectory cannot miss the inner hyperbolic surface
          // for z > 0, if the normal is correct.
          //
          G4double dot1 = (xi*v.x() + yi*v.y())*endInnerRadius/std::sqrt(pr2);
          couldMissInner = (dot1 - halfLenZ*tanInnerStereo2*vz <= 0);
          
          if (pr2 > endInnerRadius2*(1 - 2*DBL_EPSILON) )
          {
            //
            // There is a potential leak if the inner
            // surface is a cylinder
            //
            if ( (innerStereo < DBL_MIN)
              && ((std::fabs(v.x()) > DBL_MIN) || (std::fabs(v.y()) > DBL_MIN)))
              cantMissInnerCylinder = true;
          }
        }
        else
        {
          return (sigz < fHalfTol) ? 0 : q;
        }
      }
      else
      {
        G4double dotR( xi*v.x() + yi*v.y() );
        if (dotR >= 0)
        {
          //
          // Otherwise, if we are traveling outwards, we know
          // we must miss the hyperbolic surfaces also, so
          // we need not bother checking
          //
          return kInfinity;
        }
        else
        {
          //
          // This test is sufficient to ensure that the
          // trajectory cannot miss the outer hyperbolic surface
          // for z > 0, if the normal is correct.
          //
          G4double dot1 = dotR*endOuterRadius/std::sqrt(pr2);
          couldMissOuter = (dot1 - halfLenZ*tanOuterStereo2*vz>= 0);
        }
      }
    }
  }
    
  //
  // Check intersection with outer hyperbolic surface, save
  // distance to valid intersection into "best".
  //    
  G4double best = kInfinity;
  
  G4double q[2];
  G4int n = IntersectHype( p, v, outerRadius2, tanOuterStereo2, q );
  
  if (n > 0)
  {
    //
    // Potential intersection: is p on this surface?
    //
    if (pz < halfLenZ+fHalfTol)
    {
      G4double dr2 = p.x()*p.x() + p.y()*p.y() - HypeOuterRadius2(pz);
      if (std::fabs(dr2) < kCarTolerance*endOuterRadius)
      {
        //
        // Sure, but make sure we're traveling inwards at
        // this point
        //
        if (p.x()*v.x() + p.y()*v.y() - pz*tanOuterStereo2*vz < 0)
          return 0;
      }
    }
    
    //
    // We are now certain that p is not on the tolerant surface.
    // Accept only position distance q
    //
    G4int i;
    for( i=0; i<n; i++ )
    {
      if (q[i] >= 0)
      {
        //
        // Check to make sure this intersection point is
        // on the surface, but only do so if we haven't
        // checked the endplate intersection already
        //
        G4double zi = pz + q[i]*vz;
        
        if (zi < -halfLenZ) continue;
        if (zi > +halfLenZ && couldMissOuter) continue;
        
        //
        // Check normal
        //
        G4double xi = p.x() + q[i]*v.x(),
           yi = p.y() + q[i]*v.y();
           
        if (xi*v.x() + yi*v.y() - zi*tanOuterStereo2*vz > 0) continue;

        best = q[i];
        break;
      }
    }
  }
  
  if (!InnerSurfaceExists()) return best;    
  
  //
  // Check intersection with inner hyperbolic surface
  //
  n = IntersectHype( p, v, innerRadius2, tanInnerStereo2, q );  
  if (n == 0)
  {
    if (cantMissInnerCylinder) return (sigz < fHalfTol) ? 0 : -sigz/vz;
        
    return best;
  }
  
  //
  // P on this surface?
  //
  if (pz < halfLenZ+fHalfTol)
  {
    G4double dr2 = p.x()*p.x() + p.y()*p.y() - HypeInnerRadius2(pz);
    if (std::fabs(dr2) < kCarTolerance*endInnerRadius)
    {
      //
      // Sure, but make sure we're traveling outwards at
      // this point
      //
      if (p.x()*v.x() + p.y()*v.y() - pz*tanInnerStereo2*vz > 0) return 0;
    }
  }
  
  //
  // No, so only positive q is valid. Search for a valid intersection
  // that is closer than the outer intersection (if it exists)
  //
  G4int i;
  for( i=0; i<n; i++ )
  {
    if (q[i] > best) break;
    if (q[i] >= 0)
    {
      //
      // Check to make sure this intersection point is
      // on the surface, but only do so if we haven't
      // checked the endplate intersection already
      //
      G4double zi = pz + q[i]*vz;

      if (zi < -halfLenZ) continue;
      if (zi > +halfLenZ && couldMissInner) continue;

      //
      // Check normal
      //
      G4double xi = p.x() + q[i]*v.x(),
         yi = p.y() + q[i]*v.y();

      if (xi*v.x() + yi*v.y() - zi*tanOuterStereo2*vz < 0) continue;

      best = q[i];
      break;
    }
  }
    
  //
  // Done
  //
  return best;
}
 

//
// Calculate distance to shape from outside, along perpendicular direction 
// (if one exists). May be an underestimate.
//
// There are five (r,z) regions:
//    1. a point that is beyond the endcap but within the
//       endcap radii
//    2. a point with r > outer endcap radius and with
//       a z position that is beyond the cone formed by the
//       normal of the outer hyperbolic surface at the 
//       edge at which it meets the endcap. 
//    3. a point that is outside the outer surface and not in (1 or 2)
//    4. a point that is inside the inner surface and not in (5)
//    5. a point with radius < inner endcap radius and
//       with a z position beyond the cone formed by the
//       normal of the inner hyperbolic surface at the
//       edge at which it meets the endcap.
// (regions 4 and 5 only exist if there is an inner surface)
//
G4double G4Hype::DistanceToIn(const G4ThreeVector& p) const
{
  G4double absZ(std::fabs(p.z()));
  
  //
  // Check region
  //
  G4double r2 = p.x()*p.x() + p.y()*p.y();
  G4double r = std::sqrt(r2);
  
  G4double sigz = absZ - halfLenZ;
  
  if (r < endOuterRadius)
  {
    if (sigz > -fHalfTol)
    {
      if (InnerSurfaceExists())
      {
        if (r > endInnerRadius) 
          return sigz < fHalfTol ? 0 : sigz;  // Region 1
        
        G4double dr = endInnerRadius - r;
        if (sigz > dr*tanInnerStereo2)
        {
          //
          // In region 5
          //
          G4double answer = std::sqrt( dr*dr + sigz*sigz );
          return answer < fHalfTol ? 0 : answer;
        }
      }
      else
      {
        //
        // In region 1 (no inner surface)
        //
        return sigz < fHalfTol ? 0 : sigz;
      }
    }
  }
  else
  {
    G4double dr = r - endOuterRadius;
    if (sigz > -dr*tanOuterStereo2)
    {
      //
      // In region 2
      //
      G4double answer = std::sqrt( dr*dr + sigz*sigz );
      return answer < fHalfTol ? 0 : answer;
    }
  }
  
  if (InnerSurfaceExists())
  {
    if (r2 < HypeInnerRadius2(absZ)+kCarTolerance*endInnerRadius)
    {
       //
      // In region 4
      //
      G4double answer = ApproxDistInside( r,absZ,innerRadius,tanInnerStereo2 );
      return answer < fHalfTol ? 0 : answer;
    }
  }
  
  //
  // We are left by elimination with region 3
  //
  G4double answer = ApproxDistOutside( r, absZ, outerRadius, tanOuterStereo );
  return answer < fHalfTol ? 0 : answer;
}


//
// Calculate distance to surface of shape from `inside', allowing for tolerance
//
// The situation here is much simplier than DistanceToIn(p,v). For
// example, there is no need to even check whether an intersection
// point is inside the boundary of a surface, as long as all surfaces 
// are checked and the smallest distance is used.
//
G4double G4Hype::DistanceToOut( const G4ThreeVector& p, const G4ThreeVector& v,
                                const G4bool calcNorm,
                                G4bool *validNorm, G4ThreeVector *norm ) const
{
  static const G4ThreeVector normEnd1(0.0,0.0,+1.0);
  static const G4ThreeVector normEnd2(0.0,0.0,-1.0);
  
  //
  // Keep track of closest surface
  //
  G4double sBest;        // distance to
  const G4ThreeVector *nBest;    // normal vector
  G4bool vBest;        // whether "valid"

  //
  // Check endplate, taking advantage of symmetry.
  // Note that the endcap is the only surface which
  // has a "valid" normal, i.e. is a surface of which
  // the entire solid is behind.
  //
  G4double pz(p.z()), vz(v.z());
  if (vz < 0)
  {
    pz = -pz;
    vz = -vz;
    nBest = &normEnd2;
  }
  else
    nBest = &normEnd1;

  //
  // Possible intercept. Are we on the surface?
  //
  if (pz > halfLenZ-fHalfTol)
  {
    if (calcNorm) { *norm = *nBest; *validNorm = true; }
    return 0;
  }

  //
  // Nope. Get distance. Beware of zero vz.
  //
  sBest = (vz > DBL_MIN) ? (halfLenZ - pz)/vz : kInfinity;
  vBest = true;
  
  //
  // Check outer surface
  //
  G4double r2 = p.x()*p.x() + p.y()*p.y();
  
  G4double q[2];
  G4int n = IntersectHype( p, v, outerRadius2, tanOuterStereo2, q );
  
  G4ThreeVector norm1, norm2;

  if (n > 0)
  {
    //
    // We hit somewhere. Are we on the surface?
    //  
    G4double dr2 = r2 - HypeOuterRadius2(pz);
    if (std::fabs(dr2) < endOuterRadius*kCarTolerance)
    {
      G4ThreeVector normHere( p.x(), p.y(), -p.z()*tanOuterStereo2 );
      //
      // Sure. But are we going the right way?
      //
      if (normHere.dot(v) > 0)
      {
        if (calcNorm) { *norm = normHere.unit(); *validNorm = false; }
        return 0;
      }
    }
    
    //
    // Nope. Check closest positive intercept.
    //
    G4int i;
    for( i=0; i<n; i++ )
    {
      if (q[i] > sBest) break;
      if (q[i] > 0)
      {
        //
        // Make sure normal is correct (that this
        // solution is an outgoing solution)
        //
        G4ThreeVector pk(p+q[i]*v);
        norm1 = G4ThreeVector( pk.x(), pk.y(), -pk.z()*tanOuterStereo2 );
        if (norm1.dot(v) > 0)
        {
          sBest = q[i];
          nBest = &norm1;
          vBest = false;
          break;
        }
      }
    }
  }
  
  if (InnerSurfaceExists())
  {
    //
    // Check inner surface
    //
    n = IntersectHype( p, v, innerRadius2, tanInnerStereo2, q );
    if (n > 0)
    {
      //
      // On surface?
      //
      G4double dr2 = r2 - HypeInnerRadius2(pz);
      if (std::fabs(dr2) < endInnerRadius*kCarTolerance)
      {
        G4ThreeVector normHere( -p.x(), -p.y(), p.z()*tanInnerStereo2 );
        if (normHere.dot(v) > 0)
        {
          if (calcNorm)
          {
            *norm = normHere.unit();
            *validNorm = false;
          }
          return 0;
        }
      }
      
      //
      // Check closest positive
      //
      G4int i;
      for( i=0; i<n; i++ )
      {
        if (q[i] > sBest) break;
        if (q[i] > 0)
        {
          G4ThreeVector pk(p+q[i]*v);
          norm2 = G4ThreeVector( -pk.x(), -pk.y(), pk.z()*tanInnerStereo2 );
          if (norm2.dot(v) > 0)
          {
            sBest = q[i];
            nBest = &norm2;
            vBest = false;
            break;
          }
        }
      }
    }
  }
  
  //
  // Done!
  //
  if (calcNorm)
  {
    *validNorm = vBest;
    
    if (nBest == &norm1 || nBest == &norm2) 
      *norm = nBest->unit();
    else
      *norm = *nBest;
  }
  
  return sBest;
}


//
// Calculate distance (<=actual) to closest surface of shape from inside
//
// May be an underestimate
//
G4double G4Hype::DistanceToOut(const G4ThreeVector& p) const
{
  //
  // Try each surface and remember the closest
  //
  G4double absZ(std::fabs(p.z()));
  G4double r(p.perp());
  
  G4double sBest = halfLenZ - absZ;
  
  G4double tryOuter = ApproxDistInside( r, absZ, outerRadius, tanOuterStereo2 );
  if (tryOuter < sBest)
    sBest = tryOuter;
  
  if (InnerSurfaceExists())
  {
    G4double tryInner = ApproxDistOutside( r,absZ,innerRadius,tanInnerStereo );
    if (tryInner < sBest) sBest = tryInner;
  }
  
  return sBest < 0.5*kCarTolerance ? 0 : sBest;
}


//
// IntersectHype (static)
//
// Decide if and where a line intersects with a hyperbolic
// surface (of infinite extent)
//
// Arguments:
//     p       - (in) Point on trajectory
//     v       - (in) Vector along trajectory
//     r2      - (in) Square of radius at z = 0
//     tan2phi - (in) std::tan(phi)**2
//     q       - (out) Up to two points of intersection, where the
//                     intersection point is p + q*v, and if there are
//                     two intersections, q[0] < q[1]. May be negative.
// Returns:
//     The number of intersections. If 0, the trajectory misses. 
//
//
// Equation of a line:
//
//       x = x0 + q*tx      y = y0 + q*ty      z = z0 + q*tz
//
// Equation of a hyperbolic surface:
//
//       x**2 + y**2 = r**2 + (z*tanPhi)**2
//
// Solution is quadratic:
//
//  a*q**2 + b*q + c = 0
//
// where:
//
//  a = tx**2 + ty**2 - (tz*tanPhi)**2
//
//  b = 2*( x0*tx + y0*ty - z0*tz*tanPhi**2 )
//
//  c = x0**2 + y0**2 - r**2 - (z0*tanPhi)**2
//
// 
G4int G4Hype::IntersectHype( const G4ThreeVector &p, const G4ThreeVector &v, 
                             G4double r2, G4double tan2Phi, G4double ss[2] )
{
  G4double x0 = p.x(), y0 = p.y(), z0 = p.z();
  G4double tx = v.x(), ty = v.y(), tz = v.z();

  G4double a = tx*tx + ty*ty - tz*tz*tan2Phi;
  G4double b = 2*( x0*tx + y0*ty - z0*tz*tan2Phi );
  G4double c = x0*x0 + y0*y0 - r2 - z0*z0*tan2Phi;
  
  if (std::fabs(a) < DBL_MIN)
  {
    //
    // The trajectory is parallel to the asympotic limit of
    // the surface: single solution
    //
    if (std::fabs(b) < DBL_MIN) return 0;
    // Unless we travel through exact center
    
    ss[0] = c/b;
    return 1;
  }
    
  
  G4double radical = b*b - 4*a*c;
  
  if (radical < -DBL_MIN) return 0;    // No solution
  
  if (radical < DBL_MIN)
  {
    //
    // Grazes surface
    //
    ss[0] = -b/a/2.0;
    return 1;
  }
  
  radical = std::sqrt(radical);
  
  G4double q = -0.5*( b + (b < 0 ? -radical : +radical) );
  G4double sa = q/a;
  G4double sb = c/q;    
  if (sa < sb) { ss[0] = sa; ss[1] = sb; } else { ss[0] = sb; ss[1] = sa; }
  return 2;
}
  
  
//
// ApproxDistOutside (static)
//
// Find the approximate distance of a point outside
// (greater radius) of a hyperbolic surface. The distance
// must be an underestimate. It will also be nice (although
// not necesary) that the estimate is always finite no
// matter how close the point is.
//
// Our hyperbola approaches the asymptotic limit at z = +/- infinity
// to the lines r = z*tanPhi. We call these lines the 
// asymptotic limit line.
//
// We need the distance of the 2d point p(r,z) to the
// hyperbola r**2 = r0**2 + (z*tanPhi)**2. Find two
// points that bracket the true normal and use the 
// distance to the line that connects these two points.
// The first such point is z=p.z. The second point is
// the z position on the asymptotic limit line that
// contains the normal on the line through the point p.
//
G4double G4Hype::ApproxDistOutside( G4double pr, G4double pz,
                                    G4double r0, G4double tanPhi )
{
  if (tanPhi < DBL_MIN) return pr-r0;

  G4double tan2Phi = tanPhi*tanPhi;

  //
  // First point
  //
  G4double z1 = pz;
  G4double r1 = std::sqrt( r0*r0 + z1*z1*tan2Phi );
  
  //
  // Second point
  //
  G4double z2 = (pr*tanPhi + pz)/(1 + tan2Phi);
  G4double r2 = std::sqrt( r0*r0 + z2*z2*tan2Phi );
  
  //
  // Line between them
  //
  G4double dr = r2-r1;
  G4double dz = z2-z1;
  
  G4double len = std::sqrt(dr*dr + dz*dz);
  if (len < DBL_MIN)
  {
    //
    // The two points are the same?? I guess we
    // must have really bracketed the normal
    //
    dr = pr-r1;
    dz = pz-z1;
    return std::sqrt( dr*dr + dz*dz );
  }
  
  //
  // Distance
  //
  return std::fabs((pr-r1)*dz - (pz-z1)*dr)/len;
}

//
// ApproxDistInside (static)
//
// Find the approximate distance of a point inside
// of a hyperbolic surface. The distance
// must be an underestimate. It will also be nice (although
// not necesary) that the estimate is always finite no
// matter how close the point is.
//
// This estimate uses the distance to a line tangent to
// the hyperbolic function. The point of tangent is chosen
// by the z position point
//
// Assumes pr and pz are positive
//
G4double G4Hype::ApproxDistInside( G4double pr, G4double pz,
                                   G4double r0, G4double tan2Phi )
{
  if (tan2Phi < DBL_MIN) return r0 - pr;

  //
  // Corresponding position and normal on hyperbolic
  //
  G4double rh = std::sqrt( r0*r0 + pz*pz*tan2Phi );
  
  G4double dr = -rh;
  G4double dz = pz*tan2Phi;
  G4double len = std::sqrt(dr*dr + dz*dz);
  
  //
  // Answer
  //
  return std::fabs((pr-rh)*dr)/len;
}


//
// GetEntityType
//
G4GeometryType G4Hype::GetEntityType() const
{
  return G4String("G4Hype");
}


//
// Clone
//
G4VSolid* G4Hype::Clone() const
{
  return new G4Hype(*this);
}


//
// GetCubicVolume
//
G4double G4Hype::GetCubicVolume()
{
  if(fCubicVolume != 0.) {;}
    else { fCubicVolume = G4VSolid::GetCubicVolume(); }
  return fCubicVolume;
}


//
// GetSurfaceArea
//
G4double G4Hype::GetSurfaceArea()
{
  if(fSurfaceArea != 0.) {;}
  else   { fSurfaceArea = G4VSolid::GetSurfaceArea(); }
  return fSurfaceArea;
}


//
// Stream object contents to an output stream
//
std::ostream& G4Hype::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Hype\n"
     << " Parameters: \n"
     << "    half length Z: " << halfLenZ/mm << " mm \n"
     << "    inner radius : " << innerRadius/mm << " mm \n"
     << "    outer radius : " << outerRadius/mm << " mm \n"
     << "    inner stereo angle : " << innerStereo/degree << " degrees \n"
     << "    outer stereo angle : " << outerStereo/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}



//
// GetPointOnSurface
//
G4ThreeVector G4Hype::GetPointOnSurface() const
{
  G4double xRand, yRand, zRand, r2 , aOne, aTwo, aThree, chose, sinhu;
  G4double phi, cosphi, sinphi, rBar2Out, rBar2In, alpha, t, rOut, rIn2, rOut2;

  // we use the formula of the area of a surface of revolution to compute 
  // the areas, using the equation of the hyperbola:
  // x^2 + y^2 = (z*tanphi)^2 + r^2

  rBar2Out = outerRadius2;
  alpha = 2.*pi*rBar2Out*std::cos(outerStereo)/tanOuterStereo;
  t     = halfLenZ*tanOuterStereo/(outerRadius*std::cos(outerStereo));
  t     = std::log(t+std::sqrt(sqr(t)+1));
  aOne  = std::fabs(2.*alpha*(std::sinh(2.*t)/4.+t/2.));


  rBar2In = innerRadius2;
  alpha = 2.*pi*rBar2In*std::cos(innerStereo)/tanInnerStereo;
  t     = halfLenZ*tanInnerStereo/(innerRadius*std::cos(innerStereo));
  t     = std::log(t+std::sqrt(sqr(t)+1));
  aTwo  = std::fabs(2.*alpha*(std::sinh(2.*t)/4.+t/2.));

  aThree = pi*((outerRadius2+sqr(halfLenZ*tanOuterStereo)
              -(innerRadius2+sqr(halfLenZ*tanInnerStereo))));
  
  if(outerStereo == 0.) {aOne = std::fabs(2.*pi*outerRadius*2.*halfLenZ);}
  if(innerStereo == 0.) {aTwo = std::fabs(2.*pi*innerRadius*2.*halfLenZ);}
  
  phi = G4RandFlat::shoot(0.,2.*pi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);
  sinhu = G4RandFlat::shoot(-1.*halfLenZ*tanOuterStereo/outerRadius,
                          halfLenZ*tanOuterStereo/outerRadius);

  chose = G4RandFlat::shoot(0.,aOne+aTwo+2.*aThree);
  if(chose>=0. && chose < aOne)
  {
    if(outerStereo != 0.)
    {
      zRand = outerRadius*sinhu/tanOuterStereo;
      xRand = std::sqrt(sqr(sinhu)+1)*outerRadius*cosphi;
      yRand = std::sqrt(sqr(sinhu)+1)*outerRadius*sinphi;
      return G4ThreeVector (xRand, yRand, zRand);
    }
    else
    {
      return G4ThreeVector(outerRadius*cosphi,outerRadius*sinphi,
                           G4RandFlat::shoot(-halfLenZ,halfLenZ));
    }
  }
  else if(chose>=aOne && chose<aOne+aTwo)
  {
    if(innerStereo != 0.)
    {
      sinhu = G4RandFlat::shoot(-1.*halfLenZ*tanInnerStereo/innerRadius,
                                halfLenZ*tanInnerStereo/innerRadius);
      zRand = innerRadius*sinhu/tanInnerStereo;
      xRand = std::sqrt(sqr(sinhu)+1)*innerRadius*cosphi;
      yRand = std::sqrt(sqr(sinhu)+1)*innerRadius*sinphi;
      return G4ThreeVector (xRand, yRand, zRand);
    }
    else 
    {
      return G4ThreeVector(innerRadius*cosphi,innerRadius*sinphi,
                           G4RandFlat::shoot(-1.*halfLenZ,halfLenZ));
    }
  }
  else if(chose>=aOne+aTwo && chose<aOne+aTwo+aThree)
  {
    rIn2  = innerRadius2+tanInnerStereo2*halfLenZ*halfLenZ;
    rOut2 = outerRadius2+tanOuterStereo2*halfLenZ*halfLenZ;
    rOut  = std::sqrt(rOut2) ;
 
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      xRand = G4RandFlat::shoot(-rOut,rOut) ;
      yRand = G4RandFlat::shoot(-rOut,rOut) ;
      r2 = xRand*xRand + yRand*yRand ;
    } while ( ! ( r2 >= rIn2 && r2 <= rOut2 ) ) ;

    zRand = halfLenZ;
    return G4ThreeVector (xRand, yRand, zRand);
  }
  else
  {
    rIn2  = innerRadius2+tanInnerStereo2*halfLenZ*halfLenZ;
    rOut2 = outerRadius2+tanOuterStereo2*halfLenZ*halfLenZ;
    rOut  = std::sqrt(rOut2) ;
 
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      xRand = G4RandFlat::shoot(-rOut,rOut) ;
      yRand = G4RandFlat::shoot(-rOut,rOut) ;
      r2 = xRand*xRand + yRand*yRand ;
    } while ( ! ( r2 >= rIn2 && r2 <= rOut2 ) ) ;

    zRand = -1.*halfLenZ;
    return G4ThreeVector (xRand, yRand, zRand);
  }
}


//
// DescribeYourselfTo
//
void G4Hype::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}


//
// GetExtent
//
G4VisExtent G4Hype::GetExtent() const 
{
  // Define the sides of the box into which the G4Tubs instance would fit.
  //
  return G4VisExtent( -endOuterRadius, endOuterRadius, 
                      -endOuterRadius, endOuterRadius, 
                      -halfLenZ, halfLenZ );
}


//
// CreatePolyhedron
//
G4Polyhedron* G4Hype::CreatePolyhedron() const 
{
   return new G4PolyhedronHype(innerRadius, outerRadius,
                               tanInnerStereo2, tanOuterStereo2, halfLenZ);
}


//
// GetPolyhedron
//
G4Polyhedron* G4Hype::GetPolyhedron () const
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


//
//  asinh
//
G4double G4Hype::asinh(G4double arg)
{
  return std::log(arg+std::sqrt(sqr(arg)+1));
}
