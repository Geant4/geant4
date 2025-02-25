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
// G4GenericTrap implementation
//
// Authors:
//   Tatiana Nikitina, CERN; Ivana Hrivnacova, IPN Orsay
//   Adapted from Root Arb8 implementation by Andrei Gheata, CERN
//
//   27.05.2024 - Evgueni Tcherniaev, complete revision, speed up
// --------------------------------------------------------------------

#include "G4GenericTrap.hh"

#if !defined(G4GEOM_USE_UGENERICTRAP)

#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4QuickRand.hh"

#include "G4GeomTools.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex  = G4MUTEX_INITIALIZER;
  G4Mutex lateralareaMutex = G4MUTEX_INITIALIZER;
}

////////////////////////////////////////////////////////////////////////
//
// Constructor
//
G4GenericTrap::G4GenericTrap(const G4String& name, G4double halfZ,
                             const std::vector<G4TwoVector>& vertices)
  : G4VSolid(name)
{
  halfTolerance = 0.5*kCarTolerance;
  CheckParameters(halfZ, vertices);
  ComputeLateralSurfaces();
  ComputeBoundingBox();
  ComputeScratchLength();
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
G4GenericTrap::G4GenericTrap(__void__& a)
  : G4VSolid(a)
{
}

////////////////////////////////////////////////////////////////////////
//
// Destructor

G4GenericTrap::~G4GenericTrap()
{
}

////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4GenericTrap::G4GenericTrap(const G4GenericTrap& rhs)
  : G4VSolid(rhs),
    halfTolerance(rhs.halfTolerance), fScratch(rhs.fScratch),
    fDz(rhs.fDz), fVertices(rhs.fVertices), fIsTwisted(rhs.fIsTwisted),
    fMinBBox(rhs.fMinBBox), fMaxBBox(rhs.fMaxBBox),
    fVisSubdivisions(rhs.fVisSubdivisions),
    fSurfaceArea(rhs.fSurfaceArea), fCubicVolume(rhs.fCubicVolume)
{
  for (auto i = 0; i < 5; ++i) { fTwist[i] = rhs.fTwist[i]; }
  for (auto i = 0; i < 4; ++i) { fDelta[i] = rhs.fDelta[i]; }
  for (auto i = 0; i < 8; ++i) { fPlane[i] = rhs.fPlane[i]; }
  for (auto i = 0; i < 4; ++i) { fSurf[i] = rhs.fSurf[i]; }
  for (auto i = 0; i < 4; ++i) { fArea[i] = rhs.fArea[i]; }
}

////////////////////////////////////////////////////////////////////////
//
// Assignment
//
G4GenericTrap& G4GenericTrap::operator=(const G4GenericTrap& rhs)
{
  // Check assignment to self
  if (this == &rhs)  { return *this; }

  // Copy base class data
  G4VSolid::operator=(rhs);

  // Copy data
  halfTolerance = rhs.halfTolerance; fScratch = rhs.fScratch;
  fDz = rhs.fDz; fVertices = rhs.fVertices; fIsTwisted = rhs.fIsTwisted;
  fMinBBox = rhs.fMinBBox; fMaxBBox = rhs.fMaxBBox;
  fVisSubdivisions = rhs.fVisSubdivisions;
  fSurfaceArea = rhs.fSurfaceArea; fCubicVolume = rhs.fCubicVolume;

  for (auto i = 0; i < 8; ++i) { fVertices[i] = rhs.fVertices[i]; }
  for (auto i = 0; i < 5; ++i) { fTwist[i] = rhs.fTwist[i]; }
  for (auto i = 0; i < 4; ++i) { fDelta[i] = rhs.fDelta[i]; }
  for (auto i = 0; i < 8; ++i) { fPlane[i] = rhs.fPlane[i]; }
  for (auto i = 0; i < 4; ++i) { fSurf[i] = rhs.fSurf[i]; }
  for (auto i = 0; i < 4; ++i) { fArea[i] = rhs.fArea[i]; }

  fRebuildPolyhedron = false;
  delete fpPolyhedron; fpPolyhedron = nullptr;

  return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Returns position of the point (inside/outside/surface)
//
EInside G4GenericTrap::Inside(const G4ThreeVector& p) const
{
  G4double px = p.x(), py = p.y(), pz = p.z();

  // intersect edges by z = pz plane
  G4TwoVector pp[4];
  G4double z = (pz + fDz);
  for (auto i = 0; i < 4; ++i) { pp[i] = fVertices[i] + fDelta[i]*z; }

  // estimate distance to the solid
  G4double dist = std::abs(pz) - fDz;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.)
    {
      G4double dd = fSurf[i].D*px + fSurf[i].E*py + fSurf[i].F*pz + fSurf[i].G;
      dist = std::max(dd, dist);
    }
    else
    {
      auto j = (i + 1)%4;
      G4TwoVector a = pp[i];
      G4TwoVector b = pp[j];
      G4double dx = b.x() - a.x();
      G4double dy = b.y() - a.y();
      G4double leng = std::sqrt(dx*dx + dy*dy);
      G4double dd = (dx*(py - a.y()) - dy*(px - a.x()))/leng;
      dist = std::max(dd, dist);
    }
  }
  return (dist > halfTolerance) ? kOutside :
    ((dist > -halfTolerance) ? kSurface : kInside);
}

////////////////////////////////////////////////////////////////////////
//
// Return unit normal to surface at p
//
G4ThreeVector G4GenericTrap::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4double halfToleranceSquared = halfTolerance*halfTolerance;
  G4double px = p.x(), py = p.y(), pz = p.z();
  G4double nx = 0, ny = 0, nz = 0;

  // intersect edges by z = pz plane
  G4TwoVector pp[4];
  G4double tz = (pz + fDz);
  for (auto i = 0; i < 4; ++i) { pp[i] = fVertices[i] + fDelta[i]*tz; }

  // check bottom and top faces
  G4double dz = std::abs(pz) - fDz;
  nz = std::copysign(G4double(std::abs(dz) <= halfTolerance), pz);

  // check lateral faces
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.)
    {
      G4double dd = fSurf[i].D*px + fSurf[i].E*py + fSurf[i].F*pz + fSurf[i].G;
      if (std::abs(dd) <= halfTolerance)
      {
        nx += fSurf[i].D;
        ny += fSurf[i].E;
        nz += fSurf[i].F;
      }
    }
    else
    {
      auto j = (i + 1)%4;
      G4TwoVector a = pp[i];
      G4TwoVector b = pp[j];
      G4double dx = b.x() - a.x();
      G4double dy = b.y() - a.y();
      G4double ll = dx*dx + dy*dy;
      G4double dd = dx*(py - a.y()) - dy*(px - a.x());
      if (dd*dd <= halfToleranceSquared*ll)
      {
        G4double x = fSurf[i].A*pz + fSurf[i].D;
        G4double y = fSurf[i].B*pz + fSurf[i].E;
        G4double z = fSurf[i].A*px + fSurf[i].B*py + 2.*fSurf[i].C*pz + fSurf[i].F;
        G4double inv = 1./std::sqrt(x*x + y*y + z*z);
        nx += x*inv;
        ny += y*inv;
        nz += z*inv;
      }
    }
  }

  // return normal
  G4double mag2 = nx*nx + ny*ny + nz*nz;
  if (mag2 == 0.) return ApproxSurfaceNormal(p); // point is not on the surface
  G4double mag = std::sqrt(mag2);
  G4double inv = (mag == 1.) ? 1. : 1./mag;
  return { nx*inv, ny*inv, nz*inv };
}

////////////////////////////////////////////////////////////////////////
//
// Find surface nearest to the point and return corresponding normal
// Normally this method should not be called
//
G4ThreeVector G4GenericTrap::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
#ifdef G4SPECSDEBUG
  std::ostringstream message;
  message.precision(16);
  message << "Point p is not on surface of solid: " << GetName() << " !?\n"
          << "Position: " << p << " is "
          << ((Inside(p) == kInside) ? "inside" : "outside") << "\n";
  StreamInfo(message);
  G4Exception("G4GenericTrap::ApproxSurfaceNormal(p)", "GeomSolids1002",
              JustWarning, message );
#endif
  G4double px = p.x(), py = p.y(), pz = p.z();
  G4double dist = -DBL_MAX;
  auto iside = 0;

  // intersect edges by z = pz plane
  G4TwoVector pp[4];
  G4double tz = (pz + fDz);
  for (auto i = 0; i < 4; ++i) { pp[i] = fVertices[i] + fDelta[i]*tz; }

  // check lateral faces
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.)
    {
      G4double d = fSurf[i].D*px + fSurf[i].E*py + fSurf[i].F*pz + fSurf[i].G;
      if (d > dist) { dist = d; iside = i; }
    }
    else
    {
      auto j = (i + 1)%4;
      G4TwoVector a = pp[i];
      G4TwoVector b = pp[j];
      G4double dx = b.x() - a.x();
      G4double dy = b.y() - a.y();
      G4double leng = std::sqrt(dx*dx + dy*dy);
      G4double d = (dx*(py - a.y()) - dy*(px - a.x()))/leng;
      if (d > dist) { dist = d; iside = i; }
    }
  }
  // check bottom and top faces
  G4double distz = std::abs(pz) - fDz;
  if (distz >= dist) return { 0., 0., std::copysign(1., pz) };

  G4double x = fSurf[iside].A*pz + fSurf[iside].D;
  G4double y = fSurf[iside].B*pz + fSurf[iside].E;
  G4double z = fSurf[iside].A*px + fSurf[iside].B*py + 2.*fSurf[iside].C*pz + fSurf[iside].F;
  G4double inv = 1./std::sqrt(x*x + y*y + z*z);
  return { x*inv, y*inv, z*inv };
}

////////////////////////////////////////////////////////////////////////
//
// Compute distance to the surface from outside,
// return kInfinity if no intersection
//
G4double G4GenericTrap::DistanceToIn(const G4ThreeVector& p,
                                     const G4ThreeVector& v) const
{
  G4double px = p.x(), py = p.y(), pz = p.z();
  G4double vx = v.x(), vy = v.y(), vz = v.z();

  // Find intersections with the bounding polyhedron
  //
  if (std::abs(pz) - fDz >= -halfTolerance && pz*vz >= 0) { return kInfinity; }
  G4double invz = (vz == 0) ? kInfinity : -1./vz;
  G4double dz = std::copysign(fDz,invz);
  G4double xin  = (pz - dz)*invz;
  G4double xout = (pz + dz)*invz;

  // Check plane faces
  for (auto k = 0; k < 4; ++k)
  {
    if (fTwist[k] != 0) continue; // skip twisted faces
    const G4GenericTrapPlane& surf = fPlane[2*k];
    G4double cosa = surf.A*vx + surf.B*vy + surf.C*vz;
    G4double dist = surf.A*px + surf.B*py + surf.C*pz + surf.D;
    if (dist >= -halfTolerance)
    {
      if (cosa >= 0.) { return kInfinity; } // point flies away
      G4double tmp  = -dist/cosa;
      xin = std::max(tmp, xin);
    }
    else
    {
      if (cosa > 0) { xout = std::min(-dist/cosa, xout); }
      if (cosa < 0) { xin = std::max(-dist/cosa, xin); }
    }
  }

  // Check planes around twisted faces, each twisted face is bounded by two planes
  G4double tin  = xin;
  G4double tout = xout;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0) continue; // skip plane faces

    // check intersection with 1st bounding plane
    const G4GenericTrapPlane& surf1 = fPlane[2*i];
    G4double cosa = surf1.A*vx + surf1.B*vy + surf1.C*vz;
    G4double dist = surf1.A*px + surf1.B*py + surf1.C*pz + surf1.D;
    if (dist >= -halfTolerance)
    {
      if (cosa >= 0.) { return kInfinity; } // point flies away
      G4double tmp  = -dist/cosa;
      tin = std::max(tmp, tin);
    }
    else
    {
      if (cosa > 0) { tout = std::min(-dist/cosa, tout); }
      if (cosa < 0) { tin = std::max(-dist/cosa, tin); }
    }

    // check intersection with 2nd bounding plane
    const G4GenericTrapPlane& surf2 = fPlane[2*i + 1];
    cosa = surf2.A*vx + surf2.B*vy + surf2.C*vz;
    dist = surf2.A*px + surf2.B*py + surf2.C*pz + surf2.D;
    if (dist >= -halfTolerance)
    {
      if (cosa >= 0.) { return kInfinity; } // point flies away
      G4double tmp  = -dist/cosa;
      tin = std::max(tmp, tin);
    }
    else
    {
      if (cosa > 0) { tout = std::min(-dist/cosa, tout); }
      if (cosa < 0) { tin = std::max(-dist/cosa, tin); }
    }
  }
  if (tout - tin <= halfTolerance) { return kInfinity; } // touch or no hit

  // if point is outside of the bounding box
  // then move it to the surface of the bounding polyhedron
  G4double t0 = 0., x0 = px, y0 = py, z0 = pz;
  if (x0 < fMinBBox.x() - halfTolerance ||
      y0 < fMinBBox.y() - halfTolerance ||
      z0 < fMinBBox.z() - halfTolerance ||
      x0 > fMaxBBox.x() + halfTolerance ||
      y0 > fMaxBBox.y() + halfTolerance ||
      z0 > fMaxBBox.z() + halfTolerance)
  {
    t0 = tin;
    x0 += vx*t0;
    y0 += vy*t0;
    z0 += vz*t0;
   }

  // Check intersections with twisted faces
  //
  G4double ttin[2] = { DBL_MAX, DBL_MAX };
  G4double ttout[2] = { tout, tout };

  if (tin - xin < halfTolerance) ttin[0] = xin;
  if (xout - tout < halfTolerance) ttout[0] = xout;
  G4double tminimal = tin - halfTolerance;
  G4double tmaximal = tout + halfTolerance;

  constexpr G4double EPSILON = 1000.*DBL_EPSILON;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0) continue; // skip plane faces

    // twisted face, solve quadratic equation
    G4double ABC  = fSurf[i].A*vx + fSurf[i].B*vy + fSurf[i].C*vz;
    G4double ABCF = fSurf[i].A*x0 + fSurf[i].B*y0 + fSurf[i].C*z0 + fSurf[i].F;
    G4double A = ABC*vz;
    G4double B = 0.5*(fSurf[i].D*vx + fSurf[i].E*vy + ABCF*vz + ABC*z0);
    G4double C = fSurf[i].D*x0 + fSurf[i].E*y0 + ABCF*z0 + fSurf[i].G;
    if (std::abs(A) <= EPSILON)
    {
      // case 1: track is parallel to the surface
      if (std::abs(B) <= EPSILON)
      {
        // check position of the track relative to the surface,
        // for the reason of precision it's better to use (x0,y0,z0) instead of (px,py,pz)
        auto j = (i + 1)%4;
        G4double z = (z0 + fDz);
        G4TwoVector a = fVertices[i] + fDelta[i]*z;
        G4TwoVector b = fVertices[j] + fDelta[j]*z;
        G4double dx = b.x() - a.x();
        G4double dy = b.y() - a.y();
        G4double leng = std::sqrt(dx*dx + dy*dy);
        G4double dist = dx*(y0 - a.y()) - dy*(x0 - a.x());
        if (dist >= -halfTolerance*leng) { return kInfinity; }
        continue;
      }

      // case 2: single root
      G4double tmp = t0 - 0.5*C/B;
      // compute normal at the intersection point and check direction
      G4double x = px + vx*tmp;
      G4double y = py + vy*tmp;
      G4double z = pz + vz*tmp;
      const G4GenericTrapSurface& surf = fSurf[i];
      G4double nx = surf.A*z + surf.D;
      G4double ny = surf.B*z + surf.E;
      G4double nz = surf.A*x + surf.B*y + 2.*surf.C*z + surf.F;

      if (nx*vx + ny*vy + nz*vz >= 0.) // point is flying to outside
      {
        auto k = (i == 3) ? 0 : i + 1;
        G4double tz = (pz + fDz);
        G4TwoVector a = fVertices[i] + fDelta[i]*tz;
        G4TwoVector b = fVertices[k] + fDelta[k]*tz;
        G4double dx = b.x() - a.x();
        G4double dy = b.y() - a.y();
        G4double leng = std::sqrt(dx*dx + dy*dy);
        G4double dist = dx*(py - a.y()) - dy*(px - a.x());
        if (dist >= -halfTolerance*leng) { return kInfinity; } // point is flies away

        if (tmp < tminimal || tmp > tmaximal) continue;
        if (std::abs(tmp - ttout[0]) < halfTolerance) continue;
        if (tmp < ttout[0])
        {
          ttout[1] = ttout[0];
          ttout[0] = tmp;
        }
        else { ttout[1] = std::min(tmp, ttout[1]); }
      }
      else // point is flying to inside
      {
        if (tmp < tminimal || tmp > tmaximal) continue;
        if (std::abs(tmp - ttin[0]) < halfTolerance) continue;
        if (tmp < ttin[0])
        {
          ttin[1] = ttin[0];
          ttin[0] = tmp;
        }
        else { ttin[1] = std::min(tmp, ttin[1]); }
      }
      continue;
    }

    // case 3: scratch or no intersection
    G4double D = B*B - A*C;
    if (D < 0.25*fScratch*fScratch*A*A)
    {
      if (A > 0) return kInfinity;
      continue;
    }

    // case 4: two intersection points
    G4double tmp = -B - std::copysign(std::sqrt(D), B);
    G4double t1 = tmp/A + t0;
    G4double t2 = C/tmp + t0;
    G4double tsurfin = std::min(t1, t2);
    G4double tsurfout = std::max(t1, t2);
    if (A < 0) std::swap(tsurfin, tsurfout);

    if (tsurfin >= tminimal && tsurfin <= tmaximal)
    {
      if (std::abs(tsurfin - ttin[0]) >= halfTolerance)
      {
        if (tsurfin < ttin[0])
        {
          ttin[1] = ttin[0];
          ttin[0] = tsurfin;
        }
        else { ttin[1] = std::min(tsurfin, ttin[1]); }
      }
    }
    if (tsurfout >= tminimal && tsurfout <= tmaximal)
    {
      if (std::abs(tsurfout - ttout[0]) >= halfTolerance)
      {
        if (tsurfout < ttout[0])
        {
          ttout[1] = ttout[0];
          ttout[0] = tsurfout;
        }
        else { ttout[1] = std::min(tsurfout, ttout[1]); }
      }
    }
  }

  // Compute distance to In
  //
  if (ttin[0] == DBL_MAX) { return kInfinity; } // no entry point

  // single entry point
  if (ttin[1] == DBL_MAX)
  {
    G4double distin = ttin[0];
    G4double distout = (distin >= ttout[0] - halfTolerance) ? ttout[1] : ttout[0];
    G4double dist = (distout <= halfTolerance || distout - distin <= halfTolerance) ? kInfinity : distin;
    return (dist < halfTolerance) ? 0. : dist;
  }

  // two entry points
  if (ttin[1] < ttout[0])
  {
    G4double distin = ttin[1], distout = ttout[0];
    G4double dist = (distout <= halfTolerance || distout - distin <= halfTolerance) ? kInfinity : distin;
    return (dist < halfTolerance) ? 0. : dist;
  }

  // check 1st pair of in-out points
  G4double distin1 = ttin[0], distout1 = ttout[0];
  G4double dist1 = (distout1 <= halfTolerance || distout1 - distin1 <= halfTolerance) ? kInfinity : distin1;
  if (dist1 != kInfinity) { return (dist1 < halfTolerance) ? 0. : dist1; }

  // check 2nd pair of in-out points
  G4double distin2 = ttin[1], distout2 = ttout[1];
  G4double dist2 = (distout2 <= halfTolerance || distout2 - distin2 <= halfTolerance) ? kInfinity : distin2;
  return (dist2 < halfTolerance) ? 0. : dist2;
}

////////////////////////////////////////////////////////////////////////
//
// Estimate safety distance to the surface from outside
//
G4double G4GenericTrap::DistanceToIn(const G4ThreeVector& p) const
{
  G4double px = p.x(), py = p.y(), pz = p.z();

  // intersect edges by z = pz plane
  G4TwoVector pp[4];
  G4double z = (pz + fDz);
  for (auto i = 0; i < 4; ++i) { pp[i] = fVertices[i] + fDelta[i]*z; }

  // estimate distance to the solid
  G4double dist = std::abs(pz) - fDz;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.)
    {
      G4double dd = fSurf[i].D*px + fSurf[i].E*py + fSurf[i].F*pz + fSurf[i].G;
      dist = std::max(dd, dist);
    }
    else
    {
      // comptute distance to the wedge (two planes) in front of the surface
      auto j = (i + 1)%4;
      G4double cx = pp[j].x() - pp[i].x();
      G4double cy = pp[j].y() - pp[i].y();
      G4double d = (pp[i].x() - px)*cy + (py - pp[i].y())*cx;
      G4ThreeVector na(-cy, cx, fDelta[i].x()*cy - fDelta[i].y()*cx);
      G4ThreeVector nb(-cy, cx, fDelta[j].x()*cy - fDelta[j].y()*cx);
      G4double amag2 = na.mag2();
      G4double bmag2 = nb.mag2();
      G4double distab = (amag2 > bmag2) ? d/std::sqrt(amag2) : d/std::sqrt(bmag2); // d > 0
      dist = std::max(distab, dist);
    }
  }
  return (dist > 0.) ? dist : 0.; // return safety distance
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface from inside
//
G4double G4GenericTrap::DistanceToOut(const G4ThreeVector& p,
                                      const G4ThreeVector& v,
                                      const G4bool calcNorm,
                                            G4bool* validNorm,
                                            G4ThreeVector* n) const
{
  G4double px = p.x(), py = p.y(), pz = p.z();
  G4double vx = v.x(), vy = v.y(), vz = v.z();

  // Check intersections with plane faces
  //
  if ((std::abs(pz) - fDz) >= -halfTolerance && pz*vz > 0.)
  {
    if (calcNorm)
    {
      *validNorm = true;
      n->set(0., 0., std::copysign(1., pz));
    }
    return 0.;
  }
  G4double tout = (vz == 0) ? DBL_MAX : (std::copysign(fDz, vz) - pz)/vz;
  G4int iface = (vz < 0) ? -4 : -2; // little trick for z-normal: (-4+3)=-1, (-2+3)=+1

  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] != 0) continue; // skip twisted faces
    const G4GenericTrapPlane& surf = fPlane[2*i];
    G4double cosa = surf.A*vx + surf.B*vy + surf.C*vz;
    if (cosa > 0)
    {
      G4double dist = surf.A*px + surf.B*py + surf.C*pz + surf.D;
      if (dist >= -halfTolerance)
      {
        if (calcNorm)
        {
           *validNorm = true;
           n->set(surf.A, surf.B, surf.C);
        }
        return 0.;
      }
      G4double tmp = -dist/cosa;
      if (tout > tmp) { tout = tmp; iface = i; }
    }
  }

  // Check intersections with twisted faces
  //
  constexpr G4double EPSILON = 1000.*DBL_EPSILON;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0) continue; // skip plane faces

    // twisted face, solve quadratic equation
    const G4GenericTrapSurface& surf = fSurf[i];
    G4double ABC  = surf.A*vx + surf.B*vy + surf.C*vz;
    G4double ABCF = surf.A*px + surf.B*py + surf.C*pz + surf.F;
    G4double A = ABC*vz;
    G4double B = 0.5*(surf.D*vx + surf.E*vy + ABCF*vz + ABC*pz);
    G4double C = surf.D*px + surf.E*py + ABCF*pz + surf.G;

    if (std::abs(A) <= EPSILON)
    {
      // case 1: track is parallel to the surface
      if (std::abs(B) <= EPSILON) { continue; }

      // case 2: single root
      G4double tmp = -0.5*C/B;
      // compute normal at intersection point and check direction
      G4double x = px + vx*tmp;
      G4double y = py + vy*tmp;
      G4double z = pz + vz*tmp;
      G4double nx = surf.A*z + surf.D;
      G4double ny = surf.B*z + surf.E;
      G4double nz = surf.A*x + surf.B*y + 2.*surf.C*z + surf.F;

      if (nx*vx + ny*vy + nz*vz > 0.) // point is flying to outside
      {
        auto k = (i + 1)%4;
        G4double tz = (pz + fDz);
        G4TwoVector a = fVertices[i] + fDelta[i]*tz;
        G4TwoVector b = fVertices[k] + fDelta[k]*tz;
        G4double dx = b.x() - a.x();
        G4double dy = b.y() - a.y();
        G4double leng = std::sqrt(dx*dx + dy*dy);
        G4double dist = dx*(py - a.y()) - dy*(px - a.x());
        if (dist >= -halfTolerance*leng) // point is on the surface
        {
          if (calcNorm)
          {
            *validNorm = false;
            G4double normx = surf.A*pz + surf.D;
            G4double normy = surf.B*pz + surf.E;
            G4double normz = surf.A*px + surf.B*py + 2.*surf.C*pz + surf.F;
            G4double inv = 1./std::sqrt(normx*normx + normy*normy + normz*normz);
            n->set(normx*inv, normy*inv, normz*inv);
          }
          return 0.;
        }
        if (tout > tmp) { tout = tmp; iface = i; }
      }
      continue;
    }

    // case 3: scratch or no intersection
    G4double D = B*B - A*C;
    if (D < 0.25*fScratch*fScratch*A*A)
    {
      // check position of the point
      auto j = (i + 1)%4;
      G4double tz = pz + fDz;
      G4TwoVector a = fVertices[i] + fDelta[i]*tz;
      G4TwoVector b = fVertices[j] + fDelta[j]*tz;
      G4double dx = b.x() - a.x();
      G4double dy = b.y() - a.y();
      G4double leng = std::sqrt(dx*dx + dy*dy);
      G4double dist = dx*(py - a.y()) - dy*(px - a.x());
      if  (dist <= -halfTolerance*leng) { continue; } // point is inside
      if  (A > 0 || dist > halfTolerance*leng) // convex surface (or point is outside)
      {
        if (calcNorm)
        {
          *validNorm = false;
          G4double nx = surf.A*pz + surf.D;
          G4double ny = surf.B*pz + surf.E;
          G4double nz = surf.A*px + surf.B*py + 2.*surf.C*pz + surf.F;
          G4double inv = 1./std::sqrt(nx*nx + ny*ny + nz*nz);
          n->set(nx*inv, ny*inv, nz*inv);
        }
        return 0.;
      }
      continue;
    }

    // case 4: two intersection points
    G4double tmp = -B - std::copysign(std::sqrt(D), B);
    G4double t1 = tmp / A;
    G4double t2 = C / tmp;
    G4double tmin = std::min(t1, t2);
    G4double tmax = std::max(t1, t2);

    if (A < 0) // concave profile: tmin(out) -> tmax(in)
    {
      if (std::abs(tmax) < std::abs(tmin)) continue; // point flies inside
      if (tmin <= halfTolerance) // point is on external side of the surface
      {
        G4double t = 0.5*(tmin + tmax);
        G4double x = px + vx*t;
        G4double y = py + vy*t;
        G4double z = pz + vz*t;

        auto j = (i + 1)%4;
        G4double tz = z + fDz;
        G4TwoVector a = fVertices[i] + fDelta[i]*tz;
        G4TwoVector b = fVertices[j] + fDelta[j]*tz;
        G4double dx = b.x() - a.x();
        G4double dy = b.y() - a.y();
        G4double leng = std::sqrt(dx*dx + dy*dy);
        G4double dist = dx*(y - a.y()) - dy*(x - a.x());
        if  (dist <= -halfTolerance*leng) continue; // scratch
        if (calcNorm)
        {
          *validNorm = false;
          G4double nx = surf.A*pz + surf.D;
          G4double ny = surf.B*pz + surf.E;
          G4double nz = surf.A*px + surf.B*py + 2.*surf.C*pz + surf.F;
          G4double inv = 1./std::sqrt(nx*nx + ny*ny + nz*nz);
          n->set(nx*inv, ny*inv, nz*inv);
        }
        return 0.;
      }
      if (tout > tmin) { tout = tmin; iface = i; }
    }
    else // convex profile: tmin(in) -> tmax(out)
    {
      if (tmax < halfTolerance) // point is on the surface
      {
        if (calcNorm)
        {
          *validNorm = false;
          G4double nx = surf.A*pz + surf.D;
          G4double ny = surf.B*pz + surf.E;
          G4double nz = surf.A*px + surf.B*py + 2.*surf.C*pz + surf.F;
          G4double inv = 1./std::sqrt(nx*nx + ny*ny + nz*nz);
          n->set(nx*inv, ny*inv, nz*inv);
        }
        return 0.;
      }
      if (tout > tmax) { tout = tmax; iface = i; }
    }
  }

  // Compute normal, if required, and return distance to out
  //
  if (tout < halfTolerance) tout = 0.;
  if (calcNorm)
  {
    if (iface < 0)
    {
      *validNorm = true;
      n->set(0, 0, iface + 3); // little trick: (-4+3)=-1, (-2+3)=+1
    }
    else
    {
      const G4GenericTrapSurface& surf = fSurf[iface];
      if (fTwist[iface] == 0)
      {
        *validNorm = true;
        n->set(surf.D, surf.E, surf.F);
      }
      else
      {
        *validNorm = false;
        G4double x = px + vx*tout;
        G4double y = py + vy*tout;
        G4double z = pz + vz*tout;
        G4double nx = surf.A*z + surf.D;
        G4double ny = surf.B*z + surf.E;
        G4double nz = surf.A*x + surf.B*y + 2.*surf.C*z + surf.F;
        G4double inv = 1./std::sqrt(nx*nx + ny*ny + nz*nz);
        n->set(nx*inv, ny*inv, nz*inv);
      }
    }
  }
  return tout;
}

////////////////////////////////////////////////////////////////////////
//
// Estimate safety distance to the surface from inside
//
G4double G4GenericTrap::DistanceToOut(const G4ThreeVector& p) const
{
  G4double px = p.x(), py = p.y(), pz = p.z();

  // intersect edges by z = pz plane
  G4TwoVector pp[4];
  G4double z = (pz + fDz);
  for (auto i = 0; i < 4; ++i) { pp[i] = fVertices[i] + fDelta[i]*z; }

  // estimate distance to the solid
  G4double dist =  std::abs(pz) - fDz;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.)
    {
      G4double dd = fSurf[i].D*px + fSurf[i].E*py + fSurf[i].F*pz + fSurf[i].G;
      dist = std::max(dd, dist);
    }
    else
    {
      // comptute distance to the wedge (two planes) in front of the surface
      auto j = (i + 1)%4;
      G4double cx = pp[j].x() - pp[i].x();
      G4double cy = pp[j].y() - pp[i].y();
      G4double d = (pp[i].x() - px)*cy + (py - pp[i].y())*cx;
      G4ThreeVector na(-cy, cx, fDelta[i].x()*cy - fDelta[i].y()*cx);
      G4ThreeVector nb(-cy, cx, fDelta[j].x()*cy - fDelta[j].y()*cx);
      G4double amag2 = na.mag2();
      G4double bmag2 = nb.mag2();
      G4double distab = (amag2 > bmag2) ? d/std::sqrt(amag2) : d/std::sqrt(bmag2); // d < 0
      dist = std::max(distab, dist);
    }
  }
  return (dist < 0.) ? -dist : 0.; // return safety distance
}

////////////////////////////////////////////////////////////////////////
//
// Return bounding box
//
void G4GenericTrap::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  pMin = fMinBBox;
  pMax = fMaxBBox;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limits
//
G4bool
G4GenericTrap::CalculateExtent(const EAxis pAxis,
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
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = pMin < pMax;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  // To build the bounding envelope with plane faces, each lateral face of
  // the trapezoid is subdivided in two triangles. Subdivision is done by
  // duplication of vertices in the bases in a way that the envelope be
  // a convex polyhedron (some faces of the envelope can be degenerated)
  //
  G4double dz = GetZHalfLength();
  G4ThreeVectorList baseA(8), baseB(8);
  for (G4int i = 0; i < 4; ++i)
  {
    G4TwoVector va = GetVertex(i);
    G4TwoVector vb = GetVertex(i+4);
    baseA[2*i].set(va.x(), va.y(),-dz);
    baseB[2*i].set(vb.x(), vb.y(), dz);
  }
  for (G4int i = 0; i < 4; ++i)
  {
    G4int k1 = 2*i, k2 = (2*i + 2)%8;
    G4double ax = (baseA[k2].x() - baseA[k1].x());
    G4double ay = (baseA[k2].y() - baseA[k1].y());
    G4double bx = (baseB[k2].x() - baseB[k1].x());
    G4double by = (baseB[k2].y() - baseB[k1].y());
    G4double znorm = ax*by - ay*bx;
    baseA[k1+1] = (znorm < 0.0) ? baseA[k2] : baseA[k1];
    baseB[k1+1] = (znorm < 0.0) ? baseB[k1] : baseB[k2];
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin, bmax, polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

////////////////////////////////////////////////////////////////////////
//
// Return type of the solid
//
G4GeometryType G4GenericTrap::GetEntityType() const
{
  return { "G4GenericTrap" };
}

////////////////////////////////////////////////////////////////////////
//
// Return true if not twisted, false otherwise
//
G4bool G4GenericTrap::IsFaceted() const
{
  return (!fIsTwisted);
}

////////////////////////////////////////////////////////////////////////
//
// Make a clone of the solid
//
G4VSolid* G4GenericTrap::Clone() const
{
  return new G4GenericTrap(*this);
}

////////////////////////////////////////////////////////////////////////
//
// Write parameters of the solid to the specified output stream
//
std::ostream& G4GenericTrap::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << "Solid geometry type: " << GetEntityType() << "\n"
     << "   half length Z: " << fDz/mm << "\n"
     << "   list of vertices:\n";
  for (G4int i = 0; i < 8; ++i)
  {
    os << "    #" << i << " " << fVertices[i] << "\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

////////////////////////////////////////////////////////////////////////
//
// Pick up a random point on the surface
//
G4ThreeVector G4GenericTrap::GetPointOnSurface() const
{
  if (fArea[0] + fArea[1] + fArea[2] + fArea[3] == 0.)
  {
    G4AutoLock l(&lateralareaMutex);
    fArea[0] = GetLateralFaceArea(0);
    fArea[1] = GetLateralFaceArea(1);
    fArea[2] = GetLateralFaceArea(2);
    fArea[3] = GetLateralFaceArea(3);
    l.unlock();
  }

  constexpr G4int iface[6][4] =
    { {0,1,2,3}, {0,4,5,1}, {1,5,6,2}, {2,6,7,3}, {3,7,4,0}, {4,5,6,7} };

  G4bool isTwisted[6] = {false};
  for (auto i = 0; i < 4; ++i) { isTwisted[i + 1] = (fTwist[i] != 0.0); }

  G4double ssurf[6];
  G4TwoVector A = fVertices[3] - fVertices[1];
  G4TwoVector B = fVertices[2] - fVertices[0];
  G4TwoVector C = fVertices[7] - fVertices[5];
  G4TwoVector D = fVertices[6] - fVertices[4];
  ssurf[0] = (A.x()*B.y() - A.y()*B.x())*0.5; // -fDz face
  ssurf[1] = ssurf[0] + fArea[0];
  ssurf[2] = ssurf[1] + fArea[1];
  ssurf[3] = ssurf[2] + fArea[2];
  ssurf[4] = ssurf[3] + fArea[3];
  ssurf[5] = ssurf[4] + (C.x()*D.y() - C.y()*D.x())*.5; // +fDz face

  G4double select = ssurf[5]*G4QuickRand();
  G4int k;
  for (k = 0; k < 5; ++k) { if (select <= ssurf[k]) break; }

  G4int i0 = iface[k][0];
  G4int i1 = iface[k][1];
  G4int i2 = iface[k][2];
  G4int i3 = iface[k][3];
  G4ThreeVector pp[4];
  pp[0].set(fVertices[i0].x(), fVertices[i0].y(), ((k == 5) ?  fDz : -fDz));
  pp[1].set(fVertices[i1].x(), fVertices[i1].y(), ((k == 0) ? -fDz :  fDz));
  pp[2].set(fVertices[i2].x(), fVertices[i2].y(), ((k == 0) ? -fDz :  fDz));
  pp[3].set(fVertices[i3].x(), fVertices[i3].y(), ((k == 5) ?  fDz : -fDz));

  G4ThreeVector point;
  if (isTwisted[k])
  { // twisted face, rejection sampling
    G4double maxlength = std::max((pp[2] - pp[1]).mag(), (pp[3] - pp[0]).mag());
    point = (pp[0] + pp[1] + pp[2] + pp[3])*0.25;
    for (auto i = 0; i < 10000; ++i)
    {
      G4double u = G4QuickRand();
      G4ThreeVector p0 = pp[0] + (pp[1] - pp[0])*u;
      G4ThreeVector p1 = pp[3] + (pp[2] - pp[3])*u;
      G4double v = G4QuickRand()*(maxlength/(p1 - p0).mag());
      if (v >= 1.) continue;
      point = p0 + (p1 - p0)*v;
      break;
    }
  }
  else
  { // plane face
    G4double u = G4QuickRand();
    G4double v = G4QuickRand();
    if (u + v > 1.) { u = 1. - u; v = 1. - v; }
    G4double ss = (((pp[2] - pp[0]).cross(pp[3] - pp[0])).mag())*0.5;
    G4int j = (select > ssurf[k] - ss) ? 3 : 1;
    point = pp[0] + (pp[2] - pp[0])*u + (pp[j] - pp[0])*v;
  }
  return point;
}

////////////////////////////////////////////////////////////////////////
//
// Return volume of the solid
//
G4double G4GenericTrap::GetCubicVolume()
{
  if (fCubicVolume == 0.0)
  {
    // diagonals
    G4TwoVector A = fVertices[3] - fVertices[1];
    G4TwoVector B = fVertices[2] - fVertices[0];
    G4TwoVector C = fVertices[7] - fVertices[5];
    G4TwoVector D = fVertices[6] - fVertices[4];

    // kross products
    G4double AB = A.x()*B.y() - A.y()*B.x();
    G4double CD = C.x()*D.y() - C.y()*D.x();
    G4double AD = A.x()*D.y() - A.y()*D.x();
    G4double CB = C.x()*B.y() - C.y()*B.x();

    fCubicVolume = ((AB + CD)/3. + (AD + CB)/6.)*fDz;
  }
  return fCubicVolume;
}

////////////////////////////////////////////////////////////////////////
//
// Compute lateral face area
//
G4double G4GenericTrap::GetLateralFaceArea(G4int iface) const
{
  G4int i1 = iface, i2 = (i1 + 1)%4, i3 = i1 + 4, i4 = i2 + 4;

  // plane face
  if (fTwist[iface] == 0.0)
  {
    G4ThreeVector p1(fVertices[i1].x(), fVertices[i1].y(),-fDz);
    G4ThreeVector p2(fVertices[i2].x(), fVertices[i2].y(),-fDz);
    G4ThreeVector p3(fVertices[i3].x(), fVertices[i3].y(), fDz);
    G4ThreeVector p4(fVertices[i4].x(), fVertices[i4].y(), fDz);
    return ((p4 - p1).cross(p3 - p2)).mag()*0.5;
  }

  // twisted face
  constexpr G4int NSTEP = 250;
  constexpr G4double dt = 1./NSTEP;

  G4double x21 = fVertices[i2].x() - fVertices[i1].x();
  G4double y21 = fVertices[i2].y() - fVertices[i1].y();
  G4double x31 = fVertices[i3].x() - fVertices[i1].x();
  G4double y31 = fVertices[i3].y() - fVertices[i1].y();
  G4double x42 = fVertices[i4].x() - fVertices[i2].x();
  G4double y42 = fVertices[i4].y() - fVertices[i2].y();
  G4double x43 = fVertices[i4].x() - fVertices[i3].x();
  G4double y43 = fVertices[i4].y() - fVertices[i3].y();

  G4double A = x21*y43 - y21*x43;
  G4double B0 = x21*y31 - y21*x31;
  G4double B1 = x42*y31 - y42*x31;
  G4double HH = 4*fDz*fDz;
  G4double invAA = 1./(A*A);
  G4double sqrtAA = 2.*std::abs(A);
  G4double invSqrtAA = 1./sqrtAA;

  G4double area = 0.;
  for (G4int i = 0; i < NSTEP; ++i)
  {
    G4double t = (i + 0.5)*dt;
    G4double I = y21 + (y43 - y21)*t;
    G4double J = x21 + (x43 - x21)*t;
    G4double IIJJ = HH*(I*I + J*J);
    G4double B = B1*t + B0;

    G4double aa = A*A;
    G4double bb = 2.*A*B;
    G4double cc = IIJJ + B*B;

    G4double R1 = std::sqrt(aa + bb + cc);
    G4double R0 = std::sqrt(cc);
    G4double log1 = std::log(std::abs(sqrtAA*R1 + 2.*aa + bb));
    G4double log0 = std::log(std::abs(sqrtAA*R0 + bb));

    area += 0.5*R1 + 0.25*bb*invAA*(R1 - R0) + IIJJ*invSqrtAA*(log1 - log0);
  }
  return area*dt;
}

////////////////////////////////////////////////////////////////////////
//
// Return surface area of the solid
//
G4double G4GenericTrap::GetSurfaceArea()
{
  if (fSurfaceArea == 0.0)
  {
    G4TwoVector A = fVertices[3] - fVertices[1];
    G4TwoVector B = fVertices[2] - fVertices[0];
    G4TwoVector C = fVertices[7] - fVertices[5];
    G4TwoVector D = fVertices[6] - fVertices[4];
    G4double S_bot = (A.x()*B.y() - A.y()*B.x())*0.5;
    G4double S_top = (C.x()*D.y() - C.y()*D.x())*0.5;
    fArea[0] = GetLateralFaceArea(0);
    fArea[1] = GetLateralFaceArea(1);
    fArea[2] = GetLateralFaceArea(2);
    fArea[3] = GetLateralFaceArea(3);
    fSurfaceArea = S_bot + S_top + fArea[0] + fArea[1] + fArea[2] + fArea[3];
  }
  return fSurfaceArea;
}

////////////////////////////////////////////////////////////////////////
//
// Check parameters of the solid
//
void
G4GenericTrap::CheckParameters(G4double halfZ,
                               const std::vector<G4TwoVector>& vertices)
{
  // Check number of vertices
  //
  if (vertices.size() != 8)
  {
    std::ostringstream message;
    message << "Number of vertices is " << vertices.size()
            << " instead of 8 for Solid: " << GetName() << " !";
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }

  // Check dZ
  //
  if ((fDz = halfZ) < 2.*kCarTolerance)
  {
    std::ostringstream message;
    message << "Z-dimension is too small or negative (dZ = " << fDz << " mm)"
            << " for Solid: " << GetName() << " !";
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }

  // Check order of vertices
  //
  for (auto i = 0; i < 8; ++i) { fVertices.at(i) = vertices[i]; }

  // Bottom polygon area
  G4TwoVector diag1 = fVertices[2] - fVertices[0];
  G4TwoVector diag2 = fVertices[3] - fVertices[1];
  G4double ldiagbot = std::max(diag1.mag(), diag2.mag());
  G4double zbot = diag1.x()*diag2.y() - diag1.y()*diag2.x();
  if (std::abs(zbot) < ldiagbot*kCarTolerance) zbot = 0.;

  // Top polygon area
  G4TwoVector diag3 = fVertices[6] - fVertices[4];
  G4TwoVector diag4 = fVertices[7] - fVertices[5];
  G4double ldiagtop = std::max(diag3.mag(), diag4.mag());
  G4double ztop = diag3.x()*diag4.y() - diag3.y()*diag4.x();
  if (std::abs(ztop) < ldiagtop*kCarTolerance) ztop = 0.;

  if (zbot*ztop < 0.)
  {
    std::ostringstream message;
    message << "Vertices of the bottom and top polygons are defined in opposite directions !\n";
    StreamInfo(message);
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
  if ((zbot > 0.) || (ztop > 0.))
  {
    std::swap(fVertices[1], fVertices[3]);
    std::swap(fVertices[5], fVertices[7]);
    std::ostringstream message;
    message << "Vertices re-ordered in Solid: " << GetName() << " !\n"
            << "Vertices of the bottom and top polygons must be defined in a clockwise direction.";
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids1001",
                JustWarning, message);
   }

  // Check for degeneracy
  //
  G4int ndegenerated = 0;
  for (auto i = 0; i < 4; ++i)
  {
    auto j = (i + 1)%4;
    G4double lbot = (fVertices[j] - fVertices[i]).mag();
    G4double ltop = (fVertices[j + 4] - fVertices[i + 4]).mag();
    ndegenerated += (std::max(lbot, ltop) < kCarTolerance);
  }
  if (ndegenerated > 1 ||
      GetCubicVolume() < fDz*std::max(ldiagbot, ldiagtop)*kCarTolerance)
  {
    std::ostringstream message;
    message << "Degenerated solid !\n";
    StreamInfo(message);
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }

  // Check that the polygons are convex
  //
  G4bool isConvex = true;
  for (auto i = 0; i < 4; ++i)
  {
    auto j = (i + 1)%4;
    auto k = (j + 1)%4;
    G4TwoVector edge1 = fVertices[j] - fVertices[i];
    G4TwoVector edge2 = fVertices[k] - fVertices[j];
    isConvex = ((edge1.x()*edge2.y() - edge1.y()*edge2.x()) < kCarTolerance);
    if (!isConvex) break;
    G4TwoVector edge3 = fVertices[j + 4] - fVertices[i + 4];
    G4TwoVector edge4 = fVertices[k + 4] - fVertices[j + 4];
    isConvex = ((edge3.x()*edge4.y() - edge3.y()*edge4.x()) < kCarTolerance);
    if (!isConvex) break;
  }
  if (!isConvex)
  {
    std::ostringstream message;
    message << "The bottom and top faces must be convex polygons !\n";
    StreamInfo(message);
    G4Exception("G4GenericTrap::CheckParameters()", "GeomSolids0002",
                FatalException, message);
  }
}

////////////////////////////////////////////////////////////////////////
//
// Compute surface equations and twist angles of lateral faces
//
void G4GenericTrap::ComputeLateralSurfaces()
{
  for (auto i = 0; i < 4; ++i)
  {
    auto j = (i + 1)%4;
    G4ThreeVector p1(fVertices[j].x(), fVertices[j].y(), -fDz);
    G4ThreeVector p2(fVertices[i].x(), fVertices[i].y(), -fDz);
    G4ThreeVector p3(fVertices[j + 4].x(), fVertices[j + 4].y(), fDz);
    G4ThreeVector p4(fVertices[i + 4].x(), fVertices[i + 4].y(), fDz);
    G4ThreeVector ebot = p2 - p1;
    G4ThreeVector etop = p4 - p3;
    G4double lbot = ebot.mag();
    G4double ltop = etop.mag();
    G4double zcross = ebot.x()*etop.y() - ebot.y()*etop.x();
    G4double eps = kCarTolerance*std::max(lbot,ltop);
    if (std::min(lbot, ltop) < kCarTolerance || std::abs(zcross) < eps)
    { // plane surface: Dx + Ey + Fz + G = 0
      G4ThreeVector normal;
      if (std::max(lbot, ltop) < kCarTolerance) // degenerated face
      {
        auto k = (j + 1)%4;                               //      N
        auto l = (k + 1)%4;                               //   i  |  j
        G4TwoVector vl = fVertices[l] + fVertices[l + 4]; //    +---+
        G4TwoVector vi = fVertices[i] + fVertices[i + 4]; // l /     \ k
        G4TwoVector vj = fVertices[j] + fVertices[j + 4]; //  +-------+
        G4TwoVector vk = fVertices[k] + fVertices[k + 4]; //
        G4TwoVector vij = (vi - vl).unit() + (vj - vk).unit();
        G4ThreeVector epar = (p4 + p3 - p2 - p1);
        G4ThreeVector eort = epar.cross(G4ThreeVector(vij.x(), vij.y(), 0.0));
        normal = (eort.cross(epar)).unit();
      }
      else
      {
        normal = ((p4 - p1).cross(p3 - p2)).unit();
      }
      fSurf[i].D = fPlane[2*i].A = fPlane[2*i + 1].A = normal.x();
      fSurf[i].E = fPlane[2*i].B = fPlane[2*i + 1].B = normal.y();
      fSurf[i].F = fPlane[2*i].C = fPlane[2*i + 1].C = normal.z();
      fSurf[i].G = fPlane[2*i].D = fPlane[2*i + 1].D = -normal.dot((p1 + p2 + p3 + p4)/4.);
    }
    else
    { // hyperbolic paraboloid: Axz + Byz + Czz + Dx + Ey + Fz + G = 0
      fIsTwisted = true;
      G4double angle = std::acos(ebot.dot(etop)/(lbot*ltop));
      if (angle > CLHEP::halfpi)
      {
        std::ostringstream message;
        message << "Twist on " << angle/CLHEP::deg
                << " degrees, should not be more than 90 degrees !";
        StreamInfo(message);
        G4Exception("G4GenericTrap::ComputeLateralSurfaces()", "GeomSolids0002",
                    FatalException, message);
      }
      fTwist[i] = std::copysign(angle, zcross);
      // set equation of twisted surface (hyperbolic paraboloid)
      fSurf[i].A = 2.*fDz*(p4.y() - p3.y() - p2.y() + p1.y());
      fSurf[i].B =-2.*fDz*(p4.x() - p3.x() - p2.x() + p1.x());
      fSurf[i].C = ((p4.x() - p2.x())*(p3.y() - p1.y()) - (p4.y() - p2.y())*(p3.x() - p1.x()));
      fSurf[i].D = 2.*fDz*fDz*(p4.y() - p3.y() + p2.y() - p1.y());
      fSurf[i].E =-2.*fDz*fDz*(p4.x() - p3.x() + p2.x() - p1.x());
      fSurf[i].F = 2.*fDz*(p4.x()*p3.y() - p3.x()*p4.y() - p2.x()*p1.y() + p1.x()*p2.y());
      fSurf[i].G = fDz*fDz*((p4.x() + p2.x())*(p3.y() + p1.y()) - (p3.x() + p1.x())*(p4.y() + p2.y()));
      G4double magnitude =  G4ThreeVector(fSurf[i].D, fSurf[i].E, fSurf[i].F).mag();
      if (magnitude < kCarTolerance) continue;
      fSurf[i].A /= magnitude;
      fSurf[i].B /= magnitude;
      fSurf[i].C /= magnitude;
      fSurf[i].D /= magnitude;
      fSurf[i].E /= magnitude;
      fSurf[i].F /= magnitude;
      fSurf[i].G /= magnitude;
      // set planes of bounding polyhedron
      G4ThreeVector normal1, normal2;
      G4ThreeVector c1, c2;
      if (fTwist[i] < 0.)
      {
        normal1 = ((p2 - p1).cross(p4 - p1)).unit();
        normal2 = ((p3 - p4).cross(p1 - p4)).unit();
        c1 = p1;
        c2 = p4;
      }
      else
      {
        normal1 = ((p3 - p2).cross(p1 - p2)).unit();
        normal2 = ((p2 - p3).cross(p4 - p3)).unit();
        c1 = p2;
        c2 = p3;
      }
      fPlane[2*i].A = normal1.x();
      fPlane[2*i].B = normal1.y();
      fPlane[2*i].C = normal1.z();
      fPlane[2*i].D = -normal1.dot(c1);
      fPlane[2*i + 1].A = normal2.x();
      fPlane[2*i + 1].B = normal2.y();
      fPlane[2*i + 1].C = normal2.z();
      fPlane[2*i + 1].D = -normal2.dot(c2);
    }
    fDelta[i] = (fVertices[i + 4] - fVertices[i])/(2*fDz);
  }
}

////////////////////////////////////////////////////////////////////////
//
// Set bounding box
//
void G4GenericTrap::ComputeBoundingBox()
{
  G4double minX, maxX, minY, maxY;
  minX = maxX = fVertices[0].x();
  minY = maxY = fVertices[0].y();
  for (auto i = 1; i < 8; ++i)
  {
    minX = std::min(minX, fVertices[i].x());
    maxX = std::max(maxX, fVertices[i].x());
    minY = std::min(minY, fVertices[i].y());
    maxY = std::max(maxY, fVertices[i].y());
  }
  fMinBBox = G4ThreeVector(minX, minY,-fDz);
  fMaxBBox = G4ThreeVector(maxX, maxY, fDz);

  // Check correctness of the bounding box
  //
  if (minX >= maxX || minY >= maxY || -fDz >= fDz)
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << fMinBBox
            << "\npMax = " << fMaxBBox;
    G4Exception("G4GenericTrap::ComputeBoundingBox()", "GeomSolids1001",
                JustWarning, message);
    DumpInfo();
  }
}

////////////////////////////////////////////////////////////////////////
//
// Set max length of a scratch
//
void G4GenericTrap::ComputeScratchLength()
{
  G4double scratch = kInfinity;
  for (auto i = 0; i < 4; ++i)
  {
    if (fTwist[i] == 0.) continue; // skip plane face
    auto k = (i + 1)%4;
    G4ThreeVector p1(fVertices[i].x(), fVertices[i].y(), -fDz);
    G4ThreeVector p2(fVertices[k].x(), fVertices[k].y(), -fDz);
    G4ThreeVector p3(fVertices[i + 4].x(), fVertices[i + 4].y(), fDz);
    G4ThreeVector p4(fVertices[k + 4].x(), fVertices[k + 4].y(), fDz);
    G4ThreeVector p0 = (p1 + p2 + p3 + p4)*0.25; // center of the face
    G4ThreeVector norm = SurfaceNormal(p0);
    G4ThreeVector pp[2]; // points inside and outside the surface
    pp[0] = p0 - norm * halfTolerance;
    pp[1] = p0 + norm * halfTolerance;
    G4ThreeVector vv[2]; // unit vectors along the diagonals
    vv[0] = (p4 - p1).unit();
    vv[1] = (p3 - p2).unit();
    // find intersection points and compute the scratch
    for (auto ip = 0; ip < 2; ++ip)
    {
      G4double px = pp[ip].x();
      G4double py = pp[ip].y();
      G4double pz = pp[ip].z();
      for (auto iv = 0; iv < 2; ++iv)
      {
        G4double vx = vv[iv].x();
        G4double vy = vv[iv].y();
        G4double vz = vv[iv].z();
        // solve quadratic equation
        G4double ABC  = fSurf[i].A*vx + fSurf[i].B*vy + fSurf[i].C*vz;
        G4double ABCF = fSurf[i].A*px + fSurf[i].B*py + fSurf[i].C*pz + fSurf[i].F;
        G4double A = ABC*vz;
        G4double B = 0.5*(fSurf[i].D*vx + fSurf[i].E*vy + ABCF*vz + ABC*pz);
        G4double C = fSurf[i].D*px + fSurf[i].E*py + ABCF*pz + fSurf[i].G;
        G4double D = B*B - A*C;
        if (D < 0) continue;
        G4double leng = 2.*std::sqrt(D)/std::abs(A);
        scratch = std::min(leng, scratch);
      }
    }
  }
  fScratch = std::max(kCarTolerance, scratch);
}

////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4GenericTrap::GetPolyhedron () const
{
  if ( (fpPolyhedron == nullptr)
    || fRebuildPolyhedron
    || (fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
        fpPolyhedron->GetNumberOfRotationSteps()) )
  {
    G4AutoLock l(&polyhedronMutex);
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fpPolyhedron;
}

////////////////////////////////////////////////////////////////////////
//
// Method for visualisation
//
void G4GenericTrap::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

////////////////////////////////////////////////////////////////////////
//
// Return VisExtent
//
G4VisExtent G4GenericTrap::GetExtent() const
{
  return { fMinBBox.x(), fMaxBBox.x(),
           fMinBBox.y(), fMaxBBox.y(),
           fMinBBox.z(), fMaxBBox.z() };
}

// --------------------------------------------------------------------

G4Polyhedron* G4GenericTrap::CreatePolyhedron() const
{
  // Approximation of Twisted Side
  // Construct extra Points, if Twisted Side
  //
  G4Polyhedron* polyhedron;
  G4int nVertices, nFacets;

  G4int subdivisions = 0;
  if (fIsTwisted)
  {
    if (GetVisSubdivisions() != 0)
    {
      subdivisions = GetVisSubdivisions();
    }
    else
    {
      // Estimation of Number of Subdivisions for smooth visualisation
      //
      G4double maxTwist = 0.;
      for(G4int i = 0; i < 4; ++i)
      {
        if (GetTwistAngle(i) > maxTwist) { maxTwist = GetTwistAngle(i); }
      }

      // Computes bounding vectors for the shape
      //
      G4double Dx, Dy;
      Dx = 0.5*(fMaxBBox.x() - fMinBBox.y());
      Dy = 0.5*(fMaxBBox.y() - fMinBBox.y());
      if (Dy > Dx) { Dx = Dy; }

      subdivisions = 8*G4int(maxTwist/(Dx*Dx*Dx)*fDz);
      if (subdivisions < 4)  { subdivisions = 4; }
      if (subdivisions > 30) { subdivisions = 30; }
    }
  }
  G4int sub4 = 4*subdivisions;
  nVertices = 8 + subdivisions*4;
  nFacets = 6 + subdivisions*4;
  G4double cf = 1./(subdivisions + 1);
  polyhedron = new G4Polyhedron(nVertices, nFacets);

  // Set vertices
  //
  G4int icur = 0;
  for (G4int i = 0; i < 4; ++i)
  {
    G4ThreeVector v(fVertices[i].x(),fVertices[i].y(),-fDz);
    polyhedron->SetVertex(++icur, v);
  }
  for (G4int i = 0; i < subdivisions; ++i)
  {
    for (G4int j = 0; j < 4; ++j)
    {
      G4TwoVector u = fVertices[j]+cf*(i+1)*(fVertices[j+4]-fVertices[j]);
      G4ThreeVector v(u.x(),u.y(),-fDz+cf*2*fDz*(i+1));
      polyhedron->SetVertex(++icur, v);
    }
  }
  for (G4int i = 4; i < 8; ++i)
  {
    G4ThreeVector v(fVertices[i].x(),fVertices[i].y(),fDz);
    polyhedron->SetVertex(++icur, v);
  }

  // Set facets
  //
  icur = 0;
  polyhedron->SetFacet(++icur, 1, 4, 3, 2); // Z-plane
  for (G4int i = 0; i < subdivisions + 1; ++i)
  {
    G4int is = i*4;
    polyhedron->SetFacet(++icur, 5+is, 8+is, 4+is, 1+is);
    polyhedron->SetFacet(++icur, 8+is, 7+is, 3+is, 4+is);
    polyhedron->SetFacet(++icur, 7+is, 6+is, 2+is, 3+is);
    polyhedron->SetFacet(++icur, 6+is, 5+is, 1+is, 2+is);
  }
  polyhedron->SetFacet(++icur, 5+sub4, 6+sub4, 7+sub4, 8+sub4); // Z-plane

  polyhedron->SetReferences();
  polyhedron->InvertFacets();

  return polyhedron;
}

////////////////////////////////////////////////////////////////////////
//
// Print out a warning if A has an unexpected sign
//
void G4GenericTrap::WarningSignA(const G4String& method, const G4String& icase, G4double A,
                                 const G4ThreeVector& p, const G4ThreeVector& v) const
{
  std::ostringstream message;
  message.precision(16);
  message << icase << " in " << GetName() << "\n"
          << "   p" << p << "\n"
          << "   v" << v << "\n"
          << "   A = " << A << " (is "
          << ((A < 0) ? "negative, instead of positive)" : "positive, instead of negative)") << " !?\n";
  StreamInfo(message);
  const G4String function = "G4GenericTrap::DistanceTo" + method + "(p,v)";
  G4Exception(function, "GeomSolids1002", JustWarning, message );
}

////////////////////////////////////////////////////////////////////////
//
// Print out a warning if B has an unexpected sign
//
void G4GenericTrap::WarningSignB(const G4String& method, const G4String& icase,
                                 G4double f, G4double B,
                                 const G4ThreeVector& p, const G4ThreeVector& v) const
{
  std::ostringstream message;
  message.precision(16);
  message << icase << " in " << GetName() << "\n"
          << "   p" << p << "\n"
          << "   v" << v << "\n"
          << "   f = " << f << " B = " << B << " (is "
          << ((B < 0) ? "negative, instead of positive)" : "positive, instead of negative)") << " !?\n";
  StreamInfo(message);
  const G4String function = "G4GenericTrap::DistanceTo" + method + "(p,v)";
  G4Exception(function, "GeomSolids1002", JustWarning, message );
}

////////////////////////////////////////////////////////////////////////
//
// Print out a warning in DistanceToIn(p,v)
//
void G4GenericTrap::WarningDistanceToIn(G4int k, const G4ThreeVector& p, const G4ThreeVector& v,
                                        G4double tmin, G4double tmax,
                                        const G4double ttin[2], const G4double ttout[2]) const
{
  G4String check = "";
  if (ttin[1] != DBL_MAX)
  {
    G4double tcheck = 0.5*(ttout[0] + ttin[1]);
    if (Inside(p + v*tcheck) != kOutside) check = "\n   !!! check point 0.5*(ttout[0] + ttin[1]) is NOT outside !!!";
  }

  auto position = Inside(p);
  std::ostringstream message;
  message.precision(16);
  message << k << "_Unexpected sequence of intersections in solid: " << GetName() << " !?\n"
          << "   position = " << ((position == kInside) ? "kInside" : ((position == kOutside) ? "kOutside" : "kSurface")) << "\n"
          << "   p" << p << "\n"
          << "   v" << v << "\n"
          << "   range    : [" << tmin << ", " << tmax << "]\n"
          << "   ttin[2]  : "
          << ((ttin[0] == DBL_MAX) ? kInfinity : ttin[0]) << ", "
          << ((ttin[1] == DBL_MAX) ? kInfinity : ttin[1]) << "\n"
          << "   ttout[2] : "
          << ((ttout[0] == DBL_MAX) ? kInfinity : ttout[0]) << ", "
          << ((ttout[1] == DBL_MAX) ? kInfinity : ttout[1]) << check << "\n";
  StreamInfo(message);
  G4Exception("G4GenericTrap::DistanceToIn(p,v)", "GeomSolids1002",
              JustWarning, message );
}

////////////////////////////////////////////////////////////////////////
//
// Print out a warning in DistanceToOut(p,v)
//
void G4GenericTrap::WarningDistanceToOut(const G4ThreeVector& p,
                                         const G4ThreeVector& v,
                                         G4double tout) const
{
  auto position = Inside(p);
  std::ostringstream message;
  message.precision(16);
  message << "Unexpected final tout = " << tout << " in solid: " << GetName() << " !?\n"
          << "   position = " << ((position == kInside) ? "kInside" : ((position == kOutside) ? "kOutside" : "kSurface")) << "\n"
          << "   p" << p << "\n"
          << "   v" << v << "\n";
  StreamInfo(message);
  G4Exception("G4GenericTrap::DistanceToOut(p,v)", "GeomSolids1002",
              JustWarning, message );
}

#endif
