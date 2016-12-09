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
// $Id: G4Box.cc 99469 2016-09-22 15:04:36Z gcosmo $
//
// 
//
// Implementation for G4Box class
//
//  24.06.98 - V.Grichine: insideEdge in DistanceToIn(p,v)
//  20.09.98 - V.Grichine: new algorithm of DistanceToIn(p,v)
//  07.05.00 - V.Grichine: d= DistanceToIn(p,v), if d<e/2, d=0
//  09.06.00 - V.Grichine: safety in DistanceToIn(p) against Inside(p)=kOutside
//             and information before exception in DistanceToOut(p,v,...)
//  15.11.00 - D.Williams, V.Grichine: bug fixed in CalculateExtent - change
//                                     algorithm for rotated vertices
//  23.08.16 - E.Tcherniaev: use G4BoundingEnvelope for CalculateExtent()
//  20.09.16 - E.Tcherniaev: added Extent(pmin,pmax)
// --------------------------------------------------------------------

#include "G4Box.hh"

#if !defined(G4GEOM_USE_UBOX)

#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "Randomize.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor - check & set half widths

G4Box::G4Box(const G4String& pName,
                   G4double pX,
                   G4double pY,
                   G4double pZ)
  : G4CSGSolid(pName), fDx(pX), fDy(pY), fDz(pZ)
{
  delta = 0.5*kCarTolerance;
  if ( (pX < 2*kCarTolerance)
    || (pY < 2*kCarTolerance)
    || (pZ < 2*kCarTolerance) )  // limit to thickness of surfaces
  {
    std::ostringstream message;
    message << "Dimensions too small for Solid: " << GetName() << "!" << G4endl
            << "     hX, hY, hZ = " << pX << ", " << pY << ", " << pZ;
    G4Exception("G4Box::G4Box()", "GeomSolids0002", FatalException, message);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4Box::G4Box( __void__& a )
  : G4CSGSolid(a), fDx(0.), fDy(0.), fDz(0.), delta(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Box::~G4Box()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Box::G4Box(const G4Box& rhs)
  : G4CSGSolid(rhs), fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz), delta(rhs.delta)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Box& G4Box::operator = (const G4Box& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fDx = rhs.fDx;
   fDy = rhs.fDy;
   fDz = rhs.fDz;
   delta = rhs.delta;

   return *this;
}

//////////////////////////////////////////////////////////////////////////////

void G4Box::SetXHalfLength(G4double dx)
{
  if(dx > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDx = dx;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension X too small for solid: " << GetName() << "!"
            << G4endl
            << "       hX = " << dx;
    G4Exception("G4Box::SetXHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume= 0.;
  fSurfaceArea= 0.;
  fRebuildPolyhedron = true;
} 

void G4Box::SetYHalfLength(G4double dy) 
{
  if(dy > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDy = dy;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Y too small for solid: " << GetName() << "!"
            << G4endl
            << "       hY = " << dy;
    G4Exception("G4Box::SetYHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume= 0.;
  fSurfaceArea= 0.;
  fRebuildPolyhedron = true;
} 

void G4Box::SetZHalfLength(G4double dz) 
{
  if(dz > 2*kCarTolerance)  // limit to thickness of surfaces
  {
    fDz = dz;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Z too small for solid: " << GetName() << "!"
            << G4endl
            << "       hZ = " << dz;
    G4Exception("G4Box::SetZHalfLength()", "GeomSolids0002",
                FatalException, message);
  }
  fCubicVolume= 0.;
  fSurfaceArea= 0.;
  fRebuildPolyhedron = true;
} 

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Box::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Box::Extent(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin.set(-fDx,-fDy,-fDz);
  pMax.set( fDx, fDy, fDz);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: " 
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Box::Extent()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
} 

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Box::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                    G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  Extent(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
} 

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Box::Inside(const G4ThreeVector& p) const
{
  EInside in = kOutside ;
  G4ThreeVector q(std::fabs(p.x()), std::fabs(p.y()), std::fabs(p.z()));

  if ( q.x() <= (fDx - delta) )
  {
    if (q.y() <= (fDy - delta) )
    {
      if      ( q.z() <= (fDz - delta) ) { in = kInside ;  }
      else if ( q.z() <= (fDz + delta) ) { in = kSurface ; }
    }
    else if ( q.y() <= (fDy + delta) )
    {
      if ( q.z() <= (fDz + delta) ) { in = kSurface ; }
    }
  }
  else if ( q.x() <= (fDx + delta) )
  {
    if ( q.y() <= (fDy + delta) )
    {
      if ( q.z() <= (fDz + delta) ) { in = kSurface ; }
    }
  }
  return in ;
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned

G4ThreeVector G4Box::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double distx, disty, distz ;
  G4ThreeVector norm(0.,0.,0.);

  // Calculate distances as if in 1st octant

  distx = std::fabs(std::fabs(p.x()) - fDx) ;
  disty = std::fabs(std::fabs(p.y()) - fDy) ;
  distz = std::fabs(std::fabs(p.z()) - fDz) ;

  // New code for particle on surface including edges and corners with specific
  // normals

  const G4ThreeVector nX  = G4ThreeVector( 1.0, 0,0  );
  const G4ThreeVector nmX = G4ThreeVector(-1.0, 0,0  );
  const G4ThreeVector nY  = G4ThreeVector( 0, 1.0,0  );
  const G4ThreeVector nmY = G4ThreeVector( 0,-1.0,0  );
  const G4ThreeVector nZ  = G4ThreeVector( 0, 0,  1.0);
  const G4ThreeVector nmZ = G4ThreeVector( 0, 0,- 1.0);

  G4ThreeVector normX(0.,0.,0.), normY(0.,0.,0.), normZ(0.,0.,0.);
  G4ThreeVector sumnorm(0., 0., 0.);
  G4int noSurfaces=0; 

  if (distx <= delta)         // on X/mX surface and around
  {
    noSurfaces ++; 
    if ( p.x() >= 0. )  { normX= nX ; }       // on +X surface : (1,0,0)
    else                { normX= nmX; }       //                 (-1,0,0)
    sumnorm= normX; 
  }

  if (disty <= delta)    // on one of the +Y or -Y surfaces
  {
    noSurfaces ++; 
    if ( p.y() >= 0. )  { normY= nY;  }       // on +Y surface
    else                { normY= nmY; }
    sumnorm += normY; 
  }

  if (distz <= delta)    // on one of the +Z or -Z surfaces
  {
    noSurfaces ++; 
    if ( p.z() >= 0. )  { normZ= nZ;  }       // on +Z surface
    else                { normZ= nmZ; }
    sumnorm += normZ;
  }

  static const G4double invSqrt2 = 1.0 / std::sqrt(2.0); 
  static const G4double invSqrt3 = 1.0 / std::sqrt(3.0); 

  if( noSurfaces > 0 )
  { 
    if( noSurfaces == 1 )
    { 
      norm= sumnorm; 
    }
    else
    {
      // norm = sumnorm . unit(); 
      if( noSurfaces == 2 )
      { 
        // 2 surfaces -> on edge 
        norm = invSqrt2 * sumnorm; 
      }
      else
      { 
        // 3 surfaces (on corner)
        norm = invSqrt3 * sumnorm; 
      }
    }
  }
  else
  {
#ifdef G4CSGDEBUG
     G4Exception("G4Box::SurfaceNormal(p)", "Notification", JustWarning, 
                 "Point p is not on surface !?" );
#endif 
     norm = ApproxSurfaceNormal(p);
  }
  
  return norm;
}

//////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4Box::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  G4double distx, disty, distz ;
  G4ThreeVector norm(0.,0.,0.);

  // Calculate distances as if in 1st octant

  distx = std::fabs(std::fabs(p.x()) - fDx) ;
  disty = std::fabs(std::fabs(p.y()) - fDy) ;
  distz = std::fabs(std::fabs(p.z()) - fDz) ;

  if ( distx <= disty )
  {
    if ( distx <= distz )     // Closest to X
    {
      if ( p.x() < 0 ) { norm = G4ThreeVector(-1.0,0,0) ; }
      else             { norm = G4ThreeVector( 1.0,0,0) ; }
    }
    else                      // Closest to Z
    {
      if ( p.z() < 0 ) { norm = G4ThreeVector(0,0,-1.0) ; }
      else             { norm = G4ThreeVector(0,0, 1.0) ; }
    }
  }
  else
  {
    if ( disty <= distz )      // Closest to Y
    {
      if ( p.y() < 0 ) { norm = G4ThreeVector(0,-1.0,0) ; }
      else             { norm = G4ThreeVector(0, 1.0,0) ; }
    }
    else                       // Closest to Z
    {
      if ( p.z() < 0 ) { norm = G4ThreeVector(0,0,-1.0) ; }
      else             { norm = G4ThreeVector(0,0, 1.0) ; }
    }
  }
  return norm;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possibility of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.

G4double G4Box::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{
  G4double safx, safy, safz ;
  G4double smin=0.0, sminy, sminz ; // , sminx ;
  G4double smax=kInfinity, smaxy, smaxz ; // , smaxx ;  // they always > 0
  G4double stmp ;
  G4double sOut=kInfinity, sOuty=kInfinity, sOutz=kInfinity ;

  safx = std::fabs(p.x()) - fDx ;     // minimum distance to x surface of shape
  safy = std::fabs(p.y()) - fDy ;
  safz = std::fabs(p.z()) - fDz ;

  // Will we intersect?
  // If safx/y/z is >-tol/2 the point is outside/on the box's x/y/z extent.
  // If both p.x/y/z and v.x/y/z repectively are both positive/negative,
  // travel is in a direction away from the shape.

  if (    ((p.x()*v.x() >= 0.0) && (safx > -delta)) 
       || ((p.y()*v.y() >= 0.0) && (safy > -delta))
       || ((p.z()*v.z() >= 0.0) && (safz > -delta))   ) 
  {
    return kInfinity ;  // travel away or parallel within tolerance
  }

  // Compute min / max distances for x/y/z travel:
  // X Planes

  if ( v.x() )  // != 0
  {
    stmp = 1.0/std::fabs(v.x()) ;

    if (safx >= 0.0)
    {
      smin = safx*stmp ;
      smax = (fDx+std::fabs(p.x()))*stmp ;
    }
    else
    {
      if (v.x() < 0)  { sOut = (fDx + p.x())*stmp ; }
      else            { sOut = (fDx - p.x())*stmp ; }
    }
  }

  // Y Planes

  if ( v.y() )  // != 0
  {
    stmp = 1.0/std::fabs(v.y()) ;

    if (safy >= 0.0)
    {
      sminy = safy*stmp ;
      smaxy = (fDy+std::fabs(p.y()))*stmp ;

      if (sminy > smin) { smin=sminy ; }
      if (smaxy < smax) { smax=smaxy ; }

      if (smin >= (smax-delta))
      {
        return kInfinity ;  // touch XY corner
      }
    }
    else
    {
      if (v.y() < 0)  { sOuty = (fDy + p.y())*stmp ; }
      else            { sOuty = (fDy - p.y())*stmp ; }
      if( sOuty < sOut ) { sOut = sOuty ; }
    }     
  }

  // Z planes

  if ( v.z() )  // != 0
  {
    stmp = 1.0/std::fabs(v.z()) ;

    if ( safz >= 0.0 )
    {
      sminz = safz*stmp ;
      smaxz = (fDz+std::fabs(p.z()))*stmp ;

      if (sminz > smin) { smin = sminz ; }
      if (smaxz < smax) { smax = smaxz ; }

      if (smin >= (smax-delta))
      { 
        return kInfinity ;    // touch ZX or ZY corners
      }
    }
    else
    {
      if (v.z() < 0)  { sOutz = (fDz + p.z())*stmp ; }
      else            { sOutz = (fDz - p.z())*stmp ; }
      if( sOutz < sOut ) { sOut = sOutz ; }
    }
  }

  if (sOut <= (smin + delta)) // travel over edge
  {
    return kInfinity ;
  }
  if (smin < delta)  { smin = 0.0 ; }

  return smin ;
}

//////////////////////////////////////////////////////////////////////////
// 
// Appoximate distance to box.
// Returns largest perpendicular distance to the closest x/y/z sides of
// the box, which is the most fast estimation of the shortest distance to box
// - If inside return 0

G4double G4Box::DistanceToIn(const G4ThreeVector& p) const
{
  G4double safex, safey, safez, safe = 0.0 ;

  safex = std::fabs(p.x()) - fDx ;
  safey = std::fabs(p.y()) - fDy ;
  safez = std::fabs(p.z()) - fDz ;

  if (safex > safe) { safe = safex ; }
  if (safey > safe) { safe = safey ; }
  if (safez > safe) { safe = safez ; }

  return safe ;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.
// - Eliminate one side of each pair by considering direction of v
// - when leaving a surface & v.close, return 0

G4double G4Box::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm,G4ThreeVector *n) const
{
  ESide side = kUndefined ;
  G4double pdist,stmp,snxt=kInfinity;

  if (calcNorm) { *validNorm = true ; }  // All normals are valid

  if (v.x() > 0)   // X planes
  {
    pdist = fDx - p.x() ;

    if (pdist > delta)
    {
      snxt = pdist/v.x() ;
      side = kPX ;
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(1,0,0) ; }
      return        snxt = 0 ;
    }
  }
  else if (v.x() < 0)
  {
    pdist = fDx + p.x() ;

    if (pdist > delta)
    {
      snxt = -pdist/v.x() ;
      side = kMX ;
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(-1,0,0) ; }
      return        snxt = 0 ;
    }
  }

  if (v.y() > 0)   // Y planes
  {
    pdist = fDy-p.y();

    if (pdist > delta)
    {
      stmp = pdist/v.y();

      if (stmp < snxt)
      {
        snxt = stmp;
        side = kPY;
      }
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(0,1,0) ; }
      return        snxt = 0 ;
    }
  }
  else if (v.y() < 0)
  {
    pdist = fDy + p.y() ;

    if (pdist > delta)
    {
      stmp = -pdist/v.y();

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kMY;
      }
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(0,-1,0) ; }
      return        snxt = 0 ;
    }
  }

  if (v.z() > 0)        // Z planes
  {
    pdist = fDz-p.z();

    if ( pdist > delta )
    {
      stmp = pdist/v.z();

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kPZ;
      }
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(0,0,1) ; } 
      return        snxt = 0 ;
    }
  }
  else if (v.z() < 0)
  {
    pdist = fDz + p.z();

    if ( pdist > delta )
    {
      stmp = -pdist/v.z();

      if ( stmp < snxt )
      {
        snxt = stmp;
        side = kMZ;
      }
    }
    else
    {
      if (calcNorm) { *n   = G4ThreeVector(0,0,-1) ; }
      return        snxt = 0 ;
    }
  }

  if (calcNorm)
  {      
    switch (side)
    {
      case kPX:
        *n=G4ThreeVector(1,0,0);
        break;
      case kMX:
        *n=G4ThreeVector(-1,0,0);
        break;
      case kPY:
        *n=G4ThreeVector(0,1,0);
        break;
      case kMY:
        *n=G4ThreeVector(0,-1,0);
        break;
      case kPZ:
        *n=G4ThreeVector(0,0,1);
        break;
      case kMZ:
        *n=G4ThreeVector(0,0,-1);
        break;
      default:
        G4cout << G4endl;
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Box::DistanceToOut(p,v,..)", "GeomSolids1002",
                    JustWarning, message);
        break;
    }
  }
  return snxt;
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0

G4double G4Box::DistanceToOut(const G4ThreeVector& p) const
{
  G4double safx1,safx2,safy1,safy2,safz1,safz2,safe=0.0;

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4int oldprc = G4cout.precision(16) ;
     G4cout << G4endl ;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout.precision(oldprc) ;
     G4Exception("G4Box::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, "Point p is outside !?" );
  }
#endif

  safx1 = fDx - p.x() ;
  safx2 = fDx + p.x() ;
  safy1 = fDy - p.y() ;
  safy2 = fDy + p.y() ;
  safz1 = fDz - p.z() ;
  safz2 = fDz + p.z() ;  
  
  // shortest Dist to any boundary now MIN(safx1,safx2,safy1..)

  if (safx2 < safx1) { safe = safx2; }
  else               { safe = safx1; }
  if (safy1 < safe)  { safe = safy1; }
  if (safy2 < safe)  { safe = safy2; }
  if (safz1 < safe)  { safe = safz1; }
  if (safz2 < safe)  { safe = safz2; }

  if (safe < 0) { safe = 0 ; }
  return safe ;  
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Box::GetEntityType() const
{
  return G4String("G4Box");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Box::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Box\n"
     << " Parameters: \n"
     << "    half length X: " << fDx/mm << " mm \n"
     << "    half length Y: " << fDy/mm << " mm \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Return a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface

G4ThreeVector G4Box::GetPointOnSurface() const
{
  G4double px, py, pz, select, sumS;
  G4double Sxy = fDx*fDy, Sxz = fDx*fDz, Syz = fDy*fDz;

  sumS   = Sxy + Sxz + Syz;
  select = sumS*G4UniformRand();
 
  if( select < Sxy )
  {
    px = -fDx +2*fDx*G4UniformRand();
    py = -fDy +2*fDy*G4UniformRand();

    if(G4UniformRand() > 0.5) { pz =  fDz; }
    else                      { pz = -fDz; }
  }
  else if ( ( select - Sxy ) < Sxz ) 
  {
    px = -fDx +2*fDx*G4UniformRand();
    pz = -fDz +2*fDz*G4UniformRand();

    if(G4UniformRand() > 0.5) { py =  fDy; }
    else                      { py = -fDy; }
  }
  else  
  {
    py = -fDy +2*fDy*G4UniformRand();
    pz = -fDz +2*fDz*G4UniformRand();

    if(G4UniformRand() > 0.5) { px =  fDx; }
    else                      { px = -fDx; }
  } 
  return G4ThreeVector(px,py,pz);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Box::Clone() const
{
  return new G4Box(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Box::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}

G4VisExtent G4Box::GetExtent() const 
{
  return G4VisExtent (-fDx, fDx, -fDy, fDy, -fDz, fDz);
}

G4Polyhedron* G4Box::CreatePolyhedron () const 
{
  return new G4PolyhedronBox (fDx, fDy, fDz);
}
#endif
