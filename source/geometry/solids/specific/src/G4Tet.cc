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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P                        *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// $Id: G4Tet.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// class G4Tet
//
// Implementation for G4Tet class
//
// History:
//
//  20040903 - Marcus Mendenhall, created G4Tet
//  20041101 - Marcus Mendenhall, optimized constant dot products with
//             fCdotNijk values
//  20041101 - MHM removed tracking error by clipping DistanceToOut to 0
//             for surface cases
//  20041101 - MHM many speed optimizations in if statements
//  20041101 - MHM changed vdotn comparisons to 1e-12 instead of 0.0 to
//             avoid nearly-parallel problems
//  20041102 - MHM Added extra distance into solid to DistanceToIn(p,v)
//             hit testing
//  20041102 - MHM added ability to check for degeneracy without throwing
//             G4Exception
//  20041103 - MHM removed many unused variables from class
//  20040803 - Dionysios Anninos, added GetPointOnSurface() method
//  20061112 - MHM added code for G4VSolid GetSurfaceArea()
//  20100920 - Gabriele Cosmo added copy-ctor and operator=()
//  20160924 - Evgueni Tcherniaev, use G4BoundingEnvelope for CalculateExtent()
//
// --------------------------------------------------------------------

#include "G4Tet.hh"

#if !defined(G4GEOM_USE_UTET)

const char G4Tet::CVSVers[]="$Id: G4Tet.cc 104316 2017-05-24 13:04:23Z gcosmo $";

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"

#include "G4ThreeVector.hh"

#include <cmath>

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4Tets are created.
// A Tet has all of its geometrical information precomputed

G4Tet::G4Tet(const G4String& pName,
                   G4ThreeVector anchor,
                   G4ThreeVector p2,
                   G4ThreeVector p3,
                   G4ThreeVector p4, G4bool *degeneracyFlag)
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0), warningFlag(0)
{
  // fV<x><y> is vector from vertex <y> to vertex <x>
  //
  G4ThreeVector fV21=p2-anchor;
  G4ThreeVector fV31=p3-anchor;
  G4ThreeVector fV41=p4-anchor;

  // make sure this is a correctly oriented set of points for the tetrahedron
  //
  G4double signed_vol=fV21.cross(fV31).dot(fV41);
  if(signed_vol<0.0)
  {
    G4ThreeVector temp(p4);
    p4=p3;
    p3=temp;
    temp=fV41;
    fV41=fV31;
    fV31=temp; 
  }
  fCubicVolume = std::fabs(signed_vol) / 6.;

  G4ThreeVector fV24=p2-p4;
  G4ThreeVector fV43=p4-p3;
  G4ThreeVector fV32=p3-p2;

  fXMin=std::min(std::min(std::min(anchor.x(), p2.x()),p3.x()),p4.x());
  fXMax=std::max(std::max(std::max(anchor.x(), p2.x()),p3.x()),p4.x());
  fYMin=std::min(std::min(std::min(anchor.y(), p2.y()),p3.y()),p4.y());
  fYMax=std::max(std::max(std::max(anchor.y(), p2.y()),p3.y()),p4.y());
  fZMin=std::min(std::min(std::min(anchor.z(), p2.z()),p3.z()),p4.z());
  fZMax=std::max(std::max(std::max(anchor.z(), p2.z()),p3.z()),p4.z());

  fDx=(fXMax-fXMin)*0.5; fDy=(fYMax-fYMin)*0.5; fDz=(fZMax-fZMin)*0.5;

  fMiddle=G4ThreeVector(fXMax+fXMin, fYMax+fYMin, fZMax+fZMin)*0.5;
  fMaxSize=std::max(std::max(std::max((anchor-fMiddle).mag(),
                                      (p2-fMiddle).mag()),
                             (p3-fMiddle).mag()),
                    (p4-fMiddle).mag());

  G4bool degenerate=std::fabs(signed_vol) < 1e-9*fMaxSize*fMaxSize*fMaxSize;

  if(degeneracyFlag) *degeneracyFlag=degenerate;
  else if (degenerate)
  {
    G4Exception("G4Tet::G4Tet()", "GeomSolids0002", FatalException,
                "Degenerate tetrahedron not allowed.");
  }

  fTol=1e-9*(std::fabs(fXMin)+std::fabs(fXMax)+std::fabs(fYMin)
            +std::fabs(fYMax)+std::fabs(fZMin)+std::fabs(fZMax));
  //fTol=kCarTolerance;

  fAnchor=anchor;
  fP2=p2;
  fP3=p3;
  fP4=p4;

  G4ThreeVector fCenter123=(anchor+p2+p3)*(1.0/3.0); // face center
  G4ThreeVector fCenter134=(anchor+p4+p3)*(1.0/3.0);
  G4ThreeVector fCenter142=(anchor+p4+p2)*(1.0/3.0);
  G4ThreeVector fCenter234=(p2+p3+p4)*(1.0/3.0);

  // compute area of each triangular face by cross product
  // and sum for total surface area

  G4ThreeVector normal123=fV31.cross(fV21);
  G4ThreeVector normal134=fV41.cross(fV31);
  G4ThreeVector normal142=fV21.cross(fV41);
  G4ThreeVector normal234=fV32.cross(fV43);

  fSurfaceArea=(
      normal123.mag()+
      normal134.mag()+
      normal142.mag()+
      normal234.mag()
  )/2.0;

  fNormal123=normal123.unit();
  fNormal134=normal134.unit();
  fNormal142=normal142.unit();
  fNormal234=normal234.unit();

  fCdotN123=fCenter123.dot(fNormal123);
  fCdotN134=fCenter134.dot(fNormal134);
  fCdotN142=fCenter142.dot(fNormal142);
  fCdotN234=fCenter234.dot(fNormal234);
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Tet::G4Tet( __void__& a )
  : G4VSolid(a), fCubicVolume(0.), fSurfaceArea(0.),
    fRebuildPolyhedron(false), fpPolyhedron(0),
    fAnchor(0,0,0), fP2(0,0,0), fP3(0,0,0), fP4(0,0,0), fMiddle(0,0,0),
    fNormal123(0,0,0), fNormal142(0,0,0), fNormal134(0,0,0),
    fNormal234(0,0,0), warningFlag(0),
    fCdotN123(0.), fCdotN142(0.), fCdotN134(0.), fCdotN234(0.),
    fXMin(0.), fXMax(0.), fYMin(0.), fYMax(0.), fZMin(0.), fZMax(0.),
    fDx(0.), fDy(0.), fDz(0.), fTol(0.), fMaxSize(0.)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Tet::~G4Tet()
{
  delete fpPolyhedron;  fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Tet::G4Tet(const G4Tet& rhs)
  : G4VSolid(rhs),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fRebuildPolyhedron(false), fpPolyhedron(0), fAnchor(rhs.fAnchor),
    fP2(rhs.fP2), fP3(rhs.fP3), fP4(rhs.fP4), fMiddle(rhs.fMiddle),
    fNormal123(rhs.fNormal123), fNormal142(rhs.fNormal142),
    fNormal134(rhs.fNormal134), fNormal234(rhs.fNormal234),
    warningFlag(rhs.warningFlag), fCdotN123(rhs.fCdotN123),
    fCdotN142(rhs.fCdotN142), fCdotN134(rhs.fCdotN134),
    fCdotN234(rhs.fCdotN234), fXMin(rhs.fXMin), fXMax(rhs.fXMax),
    fYMin(rhs.fYMin), fYMax(rhs.fYMax), fZMin(rhs.fZMin), fZMax(rhs.fZMax),
    fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz), fTol(rhs.fTol),
    fMaxSize(rhs.fMaxSize)
{
}


///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Tet& G4Tet::operator = (const G4Tet& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fCubicVolume = rhs.fCubicVolume; fSurfaceArea = rhs.fSurfaceArea;
   fAnchor = rhs.fAnchor;
   fP2 = rhs.fP2; fP3 = rhs.fP3; fP4 = rhs.fP4; fMiddle = rhs.fMiddle;
   fNormal123 = rhs.fNormal123; fNormal142 = rhs.fNormal142;
   fNormal134 = rhs.fNormal134; fNormal234 = rhs.fNormal234;
   warningFlag = rhs.warningFlag; fCdotN123 = rhs.fCdotN123;
   fCdotN142 = rhs.fCdotN142; fCdotN134 = rhs.fCdotN134;
   fCdotN234 = rhs.fCdotN234; fXMin = rhs.fXMin; fXMax = rhs.fXMax;
   fYMin = rhs.fYMin; fYMax = rhs.fYMax; fZMin = rhs.fZMin; fZMax = rhs.fZMax;
   fDx = rhs.fDx; fDy = rhs.fDy; fDz = rhs.fDz; fTol = rhs.fTol;
   fMaxSize = rhs.fMaxSize;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = 0;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// CheckDegeneracy

G4bool G4Tet::CheckDegeneracy( G4ThreeVector anchor,
                               G4ThreeVector p2,
                               G4ThreeVector p3,
                               G4ThreeVector p4 )
{
  G4bool result;
  G4Tet *object=new G4Tet("temp",anchor,p2,p3,p4,&result);
  delete object;
  return result;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Tet::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Tet::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin.set(fXMin,fYMin,fZMin);
  pMax.set(fXMax,fYMax,fZMax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Tet::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Tet::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                    G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  std::vector<G4ThreeVector> vec = GetVertices();

  G4ThreeVectorList anchor(1);
  anchor[0].set(vec[0].x(),vec[0].y(),vec[0].z());

  G4ThreeVectorList base(3);
  base[0].set(vec[1].x(),vec[1].y(),vec[1].z());
  base[1].set(vec[2].x(),vec[2].y(),vec[2].z());
  base[2].set(vec[3].x(),vec[3].y(),vec[3].z());

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &anchor;
  polygons[1] = &base;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Tet::Inside(const G4ThreeVector& p) const
{
  G4double r123, r134, r142, r234;

  // this is written to allow if-statement truncation so the outside test
  // (where most of the world is) can fail very quickly and efficiently

  if ( (r123=p.dot(fNormal123)-fCdotN123) > fTol ||
       (r134=p.dot(fNormal134)-fCdotN134) > fTol ||
       (r142=p.dot(fNormal142)-fCdotN142) > fTol ||
       (r234=p.dot(fNormal234)-fCdotN234) > fTol )
  {
    return kOutside; // at least one is out!
  }
  else if( (r123 < -fTol)&&(r134 < -fTol)&&(r142 < -fTol)&&(r234 < -fTol) )
  {
    return kInside; // all are definitively inside
  }
  else
  {
    return kSurface; // too close to tell
  }
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned.
// This assumes that we are looking from the inside!

G4ThreeVector G4Tet::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double r123=std::fabs(p.dot(fNormal123)-fCdotN123);
  G4double r134=std::fabs(p.dot(fNormal134)-fCdotN134);
  G4double r142=std::fabs(p.dot(fNormal142)-fCdotN142);
  G4double r234=std::fabs(p.dot(fNormal234)-fCdotN234);

  const G4double delta = 0.5*kCarTolerance;
  G4ThreeVector sumnorm(0., 0., 0.);
  G4int noSurfaces=0; 

  if (r123 <= delta)         
  {
     noSurfaces ++; 
     sumnorm= fNormal123; 
  }

  if (r134 <= delta)    
  {
     noSurfaces ++; 
     sumnorm += fNormal134; 
  }
 
  if (r142 <= delta)    
  {
     noSurfaces ++; 
     sumnorm += fNormal142;
  }
  if (r234 <= delta)    
  {
     noSurfaces ++; 
     sumnorm += fNormal234;
  }
  
  if( noSurfaces > 0 )
  { 
     if( noSurfaces == 1 )
     { 
       return sumnorm; 
     }
     else
     {
       return sumnorm.unit();
     }
  }
  else // Approximative Surface Normal
  {

    if( (r123<=r134) && (r123<=r142) && (r123<=r234) ) { return fNormal123; }
    else if ( (r134<=r142) && (r134<=r234) )           { return fNormal134; }
    else if (r142 <= r234)                             { return fNormal142; }
    return fNormal234;
  }
}
///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
// All this is very unrolled, for speed.

G4double G4Tet::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{
    G4ThreeVector vu(v.unit()), hp;
    G4double vdotn, t, tmin=kInfinity;

    G4double extraDistance=10.0*fTol; // a little ways into the solid

    vdotn=-vu.dot(fNormal123);
    if(vdotn > 1e-12)
    { // this is a candidate face, since it is pointing at us
      t=(p.dot(fNormal123)-fCdotN123)/vdotn; // #  distance to intersection
      if( (t>=-fTol) && (t<tmin) )
      { // if not true, we're going away from this face or it's not close
        hp=p+vu*(t+extraDistance); // a little beyond point of intersection
        if ( ( hp.dot(fNormal134)-fCdotN134 < 0.0 ) &&
             ( hp.dot(fNormal142)-fCdotN142 < 0.0 ) &&
             ( hp.dot(fNormal234)-fCdotN234 < 0.0 ) )
        {
          tmin=t;
        }
      }
    }

    vdotn=-vu.dot(fNormal134);
    if(vdotn > 1e-12)
    { // # this is a candidate face, since it is pointing at us
      t=(p.dot(fNormal134)-fCdotN134)/vdotn; // #  distance to intersection
      if( (t>=-fTol) && (t<tmin) )
      { // if not true, we're going away from this face
        hp=p+vu*(t+extraDistance); // a little beyond point of intersection
        if ( ( hp.dot(fNormal123)-fCdotN123 < 0.0 ) && 
             ( hp.dot(fNormal142)-fCdotN142 < 0.0 ) &&
             ( hp.dot(fNormal234)-fCdotN234 < 0.0 ) )
        {
          tmin=t;
        }
      }
    }

    vdotn=-vu.dot(fNormal142);
    if(vdotn > 1e-12)
    { // # this is a candidate face, since it is pointing at us
      t=(p.dot(fNormal142)-fCdotN142)/vdotn; // #  distance to intersection
      if( (t>=-fTol) && (t<tmin) )
      { // if not true, we're going away from this face
        hp=p+vu*(t+extraDistance); // a little beyond point of intersection
        if ( ( hp.dot(fNormal123)-fCdotN123 < 0.0 ) &&
             ( hp.dot(fNormal134)-fCdotN134 < 0.0 ) &&
             ( hp.dot(fNormal234)-fCdotN234 < 0.0 ) )
        {
          tmin=t;
        }
      }
    }

    vdotn=-vu.dot(fNormal234);
    if(vdotn > 1e-12)
    { // # this is a candidate face, since it is pointing at us
      t=(p.dot(fNormal234)-fCdotN234)/vdotn; // #  distance to intersection
      if( (t>=-fTol) && (t<tmin) )
      { // if not true, we're going away from this face
        hp=p+vu*(t+extraDistance); // a little beyond point of intersection
        if ( ( hp.dot(fNormal123)-fCdotN123 < 0.0 ) &&
             ( hp.dot(fNormal134)-fCdotN134 < 0.0 ) &&
             ( hp.dot(fNormal142)-fCdotN142 < 0.0 ) )
        {
          tmin=t;
        }
      }
    }

  return std::max(0.0,tmin);
}

//////////////////////////////////////////////////////////////////////////
// 
// Approximate distance to tet.
// returns distance to sphere centered on bounding box
// - If inside return 0

G4double G4Tet::DistanceToIn(const G4ThreeVector& p) const
{
  G4double dd=(p-fMiddle).mag() - fMaxSize - fTol;
  return std::max(0.0, dd);
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.

G4double G4Tet::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm, G4ThreeVector *n) const
{
    G4ThreeVector vu(v.unit());
    G4double t1=kInfinity,t2=kInfinity,t3=kInfinity,t4=kInfinity, vdotn, tt;

    vdotn=vu.dot(fNormal123);
    if(vdotn > 1e-12)  // #we're heading towards this face, so it is a candidate
    {
      t1=(fCdotN123-p.dot(fNormal123))/vdotn; // #  distance to intersection
    }

    vdotn=vu.dot(fNormal134);
    if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
      t2=(fCdotN134-p.dot(fNormal134))/vdotn; // #  distance to intersection
    }

    vdotn=vu.dot(fNormal142);
    if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
      t3=(fCdotN142-p.dot(fNormal142))/vdotn; // #  distance to intersection
    }

    vdotn=vu.dot(fNormal234);
    if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
    {
      t4=(fCdotN234-p.dot(fNormal234))/vdotn; // #  distance to intersection
    }

    tt=std::min(std::min(std::min(t1,t2),t3),t4);

    if (warningFlag && (tt == kInfinity || tt < -fTol))
    {
      DumpInfo();
      std::ostringstream message;
      message << "No good intersection found or already outside!?" << G4endl
              << "p = " << p / mm << "mm" << G4endl
              << "v = " << v  << G4endl
              << "t1, t2, t3, t4 (mm) "
              << t1/mm << ", " << t2/mm << ", " << t3/mm << ", " << t4/mm;
      G4Exception("G4Tet::DistanceToOut(p,v,...)", "GeomSolids1002",
                  JustWarning, message);
      if(validNorm)
      {
        *validNorm=false; // flag normal as meaningless
      }
    }
    else if(calcNorm && n)
    {
      G4ThreeVector normal;
      if(tt==t1)        { normal=fNormal123; }
      else if (tt==t2)  { normal=fNormal134; }
      else if (tt==t3)  { normal=fNormal142; }
      else if (tt==t4)  { normal=fNormal234; }
      *n=normal;
      if(validNorm) { *validNorm=true; }
    }

    return std::max(tt,0.0); // avoid tt<0.0 by a tiny bit
                             // if we are right on a face
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0

G4double G4Tet::DistanceToOut(const G4ThreeVector& p) const
{
  G4double t1,t2,t3,t4;
  t1=fCdotN123-p.dot(fNormal123); //  distance to plane, positive if inside
  t2=fCdotN134-p.dot(fNormal134); //  distance to plane
  t3=fCdotN142-p.dot(fNormal142); //  distance to plane
  t4=fCdotN234-p.dot(fNormal234); //  distance to plane

  // if any one of these is negative, we are outside,
  // so return zero in that case

  G4double tmin=std::min(std::min(std::min(t1,t2),t3),t4);
  return (tmin < fTol)? 0:tmin;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Tet::GetEntityType() const
{
  return G4String("G4Tet");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4Tet::Clone() const
{
  return new G4Tet(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Tet::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
  << "    *** Dump for solid - " << GetName() << " ***\n"
  << "    ===================================================\n"
  << " Solid type: G4Tet\n"
  << " Parameters: \n"
  << "    anchor: " << fAnchor/mm << " mm \n"
  << "    p2: " << fP2/mm << " mm \n"
  << "    p3: " << fP3/mm << " mm \n"
  << "    p4: " << fP4/mm << " mm \n"
  << "    normal123: " << fNormal123 << " \n"
  << "    normal134: " << fNormal134 << " \n"
  << "    normal142: " << fNormal142 << " \n"
  << "    normal234: " << fNormal234 << " \n"
  << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


////////////////////////////////////////////////////////////////////////
//
// GetPointOnFace
//
// Auxiliary method for get point on surface

G4ThreeVector G4Tet::GetPointOnFace(G4ThreeVector p1, G4ThreeVector p2,
                                    G4ThreeVector p3, G4double& area) const
{
  G4double lambda1,lambda2;
  G4ThreeVector v, w;

  v = p3 - p1;
  w = p1 - p2;

  lambda1 = G4RandFlat::shoot(0.,1.);
  lambda2 = G4RandFlat::shoot(0.,lambda1);

  area = 0.5*(v.cross(w)).mag();

  return (p2 + lambda1*w + lambda2*v);
}

////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Tet::GetPointOnSurface() const
{
  G4double chose,aOne,aTwo,aThree,aFour;
  G4ThreeVector p1, p2, p3, p4;
  
  p1 = GetPointOnFace(fAnchor,fP2,fP3,aOne);
  p2 = GetPointOnFace(fAnchor,fP4,fP3,aTwo);
  p3 = GetPointOnFace(fAnchor,fP4,fP2,aThree);
  p4 = GetPointOnFace(fP4,fP3,fP2,aFour);
  
  chose = G4RandFlat::shoot(0.,aOne+aTwo+aThree+aFour);
  if( (chose>=0.) && (chose <aOne) ) {return p1;}
  else if( (chose>=aOne) && (chose < aOne+aTwo) ) {return p2;}
  else if( (chose>=aOne+aTwo) && (chose<aOne+aTwo+aThree) ) {return p3;}
  return p4;
}

////////////////////////////////////////////////////////////////////////
//
// GetVertices

std::vector<G4ThreeVector> G4Tet::GetVertices() const 
{
  std::vector<G4ThreeVector> vertices(4);
  vertices[0] = fAnchor;
  vertices[1] = fP2;
  vertices[2] = fP3;
  vertices[3] = fP4;

  return vertices;
}

////////////////////////////////////////////////////////////////////////
//
// GetCubicVolume

G4double G4Tet::GetCubicVolume()
{
  return fCubicVolume;
}

////////////////////////////////////////////////////////////////////////
//
// GetSurfaceArea

G4double G4Tet::GetSurfaceArea()
{
  return fSurfaceArea;
}

// Methods for visualisation

////////////////////////////////////////////////////////////////////////
//
// DescribeYourselfTo

void G4Tet::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////////////////////////
//
// GetExtent

G4VisExtent G4Tet::GetExtent() const 
{
  return G4VisExtent (fXMin, fXMax, fYMin, fYMax, fZMin, fZMax);
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron

G4Polyhedron* G4Tet::CreatePolyhedron () const 
{
  G4Polyhedron *ph=new G4Polyhedron;
  G4double xyz[4][3];
  const G4int faces[4][4]={{1,3,2,0},{1,4,3,0},{1,2,4,0},{2,3,4,0}};
  xyz[0][0]=fAnchor.x(); xyz[0][1]=fAnchor.y(); xyz[0][2]=fAnchor.z();
  xyz[1][0]=fP2.x(); xyz[1][1]=fP2.y(); xyz[1][2]=fP2.z();
  xyz[2][0]=fP3.x(); xyz[2][1]=fP3.y(); xyz[2][2]=fP3.z();
  xyz[3][0]=fP4.x(); xyz[3][1]=fP4.y(); xyz[3][2]=fP4.z();

  ph->createPolyhedron(4,4,xyz,faces);

  return ph;
}

////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron

G4Polyhedron* G4Tet::GetPolyhedron () const
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

#endif
