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
// ********************************************************************
//
//
// $Id:$
//
// 
// class G4UAdapter
//
// Class description:
//
// Utility class for adapting VecGeom solids API to Geant4 solids. 
// NOTE: Using protected inheritance since the Adapter is supposed to
// be a G4VSolid "implemented-in-terms-of" the VecGeom UnplacedVolume_t.
// The choice of protected vs private is due to the fact that we want
// to propagate functions further down in the inheritance hierarchy.

// Author:
// 17.05.17 G.Cosmo: Adapted for G4VSolid from original G4USolids bridge
//                   class and the USolidsAdapter class in VecGeom.
// ------------------------------------------------------------------------
#ifndef G4UADAPTER_HH
#define G4UADAPTER_HH

#include "G4ThreeVector.hh"
#include "G4VSolid.hh"

// Required for inline visualization adapter functions
//
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"
#include "G4BoundingEnvelope.hh"
#include "G4AutoLock.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <base/Global.h>
#include <base/Vector3D.h>

class G4VPVParameterisation;

template <class UnplacedVolume_t>
class G4UAdapter : public G4VSolid, protected UnplacedVolume_t
{
  public:

    typedef vecgeom::Vector3D<G4double> U3Vector;

    using UnplacedVolume_t::operator delete;
    using UnplacedVolume_t::operator new;
      // VecGeom volumes have special delete/new ("AlignedBase")
      // and we need to make these functions public again

    G4UAdapter(const G4String& name)
      : G4VSolid(name), fRebuildPolyhedron(false), fPolyhedron(0)
        { kHalfTolerance = 0.5*kCarTolerance; }

    template <typename... T>
    G4UAdapter(const G4String& name, const T &... params)
      : G4VSolid(name), UnplacedVolume_t(params...),
        fRebuildPolyhedron(false), fPolyhedron(0)
        { kHalfTolerance = 0.5*kCarTolerance; }

    virtual ~G4UAdapter();

    G4bool operator==(const G4UAdapter& s) const;
      // Return true only if addresses are the same.

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                               G4double& pMin, G4double& pMax) const override;
      // Calculate the minimum and maximum extent of the solid, when under the
      // specified transform, and within the specified limits. If the solid
      // is not intersected by the region, return false, else return true.

    virtual EInside Inside(const G4ThreeVector& p) const override;
      // Returns kOutside if the point at offset p is outside the shapes
      // boundaries plus Tolerance/2, kSurface if the point is <= Tolerance/2
      // from a surface, otherwise kInside.

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;
      // Returns the outwards pointing unit normal of the shape for the
      // surface closest to the point at offset p.

    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const override;
      // Return the distance along the normalised vector v to the shape,
      // from the point at offset p. If there is no intersection, return
      // kInfinity. The first intersection resulting from `leaving' a
      // surface/volume is discarded. Hence, it is tolerant of points on
      // the surface of the shape.

    virtual G4double DistanceToIn(const G4ThreeVector& p) const override;
      // Calculate the distance to the nearest surface of a shape from an
      // outside point. The distance can be an underestimate.

    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm = false,
                                   G4bool* validNorm = 0,
                                   G4ThreeVector* n = 0) const override;
      // Return the distance along the normalised vector v to the shape,
      // from a point at an offset p inside or on the surface of the shape.
      // Intersections with surfaces, when the point is < Tolerance/2 from a
      // surface must be ignored.
      // If calcNorm==true:
      //    validNorm set true if the solid lies entirely behind or on the
      //              exiting surface.
      //    n set to exiting outwards normal vector (undefined Magnitude).
      //    validNorm set to false if the solid does not lie entirely behind
      //              or on the exiting surface
      // If calcNorm==false:
      //    validNorm and n are unused.
      //
      // Must be called as solid.DistanceToOut(p,v) or by specifying all
      // the parameters.

    virtual G4double DistanceToOut(const G4ThreeVector& p) const override;
      // Calculate the distance to the nearest surface of a shape from an
      // inside point. The distance can be an underestimate.

    virtual void ComputeDimensions(G4VPVParameterisation* p,
                                   const G4int n,
                                   const G4VPhysicalVolume* pRep) override;
      // Throw exception if ComputeDimensions called from an illegal
      // derived class.

    virtual G4double GetCubicVolume() override;
      // Returns an estimation of the solid volume in internal units.
      // This method may be overloaded by derived classes to compute the
      // exact geometrical quantity for solids where this is possible,
      // or anyway to cache the computed value.
      // Note: the computed value is NOT cached.

    virtual G4double GetSurfaceArea() override;
      // Return an estimation of the solid surface area in internal units.
      // This method may be overloaded by derived classes to compute the
      // exact geometrical quantity for solids where this is possible,
      // or anyway to cache the computed value.
      // Note: the computed value is NOT cached.

    virtual G4ThreeVector GetPointOnSurface() const override;
      // Returns a random point located on the surface of the solid.

    virtual G4GeometryType GetEntityType() const override;
      // Provide identification of the class of an object.
      // (required for persistency)

    virtual G4VSolid* Clone() const override;
      // Returns a pointer of a dynamically allocated copy of the solid.
      // Returns NULL pointer with warning in case the concrete solid does not
      // implement this method. The caller has responsibility for ownership.

    virtual std::ostream& StreamInfo(std::ostream& os) const override;
      // Dumps contents of the solid to a stream.

    virtual void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
      // A "double dispatch" function which identifies the solid
      // to the graphics scene for visualization.

    virtual G4VisExtent GetExtent()  const override;
      // Provide extent (bounding box) as possible hint to the graphics view.
    virtual G4Polyhedron* CreatePolyhedron() const override;
      // Create Polyhedron used for Visualisation
    virtual G4Polyhedron* GetPolyhedron() const override;
      // Smart access function - creates on request and stores for future
      // access.  A null pointer means "not available".

  public:  // without description

    G4UAdapter(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UAdapter(const G4UAdapter& rhs);
    G4UAdapter& operator=(const G4UAdapter& rhs);
      // Copy constructor and assignment operator.

  public:  // VecGeom overridden methods

    vecgeom::Precision
    DistanceToOut(U3Vector const &position, U3Vector const &direction,
                  vecgeom::Precision stepMax = kInfinity) const override
    {
      return UnplacedVolume_t::DistanceToOut(position, direction, stepMax);
    }

    vecgeom::EnumInside
    Inside(U3Vector const &aPoint) const override
    {
      return UnplacedVolume_t::Inside(aPoint);
    }

    vecgeom::Precision
    DistanceToIn(U3Vector const &position, U3Vector const &direction,
                 const vecgeom::Precision step_max = kInfinity) const override
    {
      return UnplacedVolume_t::DistanceToIn(position, direction, step_max);
    }

    G4bool Normal(U3Vector const &aPoint, U3Vector &aNormal) const override
    {
      return UnplacedVolume_t::Normal(aPoint, aNormal);
    }

    void Extent(U3Vector &aMin, U3Vector &aMax) const override
    {
      return UnplacedVolume_t::Extent(aMin, aMax);
    }

    U3Vector SamplePointOnSurface() const override
    {
      return UnplacedVolume_t::SamplePointOnSurface();
    }

  protected:  // data

    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fPolyhedron;

    G4double kHalfTolerance;      // Cached geometrical tolerance

    using UnplacedVolume_t::DistanceToOut;
    using UnplacedVolume_t::DistanceToIn;
};

// Inline implementations

template <class UnplacedVolume_t>
G4UAdapter<UnplacedVolume_t>::G4UAdapter(__void__& a)
  : G4VSolid(a), UnplacedVolume_t(*this),
    fRebuildPolyhedron(false), fPolyhedron(0),
    kHalfTolerance(0.5*kCarTolerance)
{
}

template <class UnplacedVolume_t>
G4UAdapter<UnplacedVolume_t>::~G4UAdapter()
{
  delete fPolyhedron; fPolyhedron = 0;
}

template <class UnplacedVolume_t>
G4bool G4UAdapter<UnplacedVolume_t>::
operator==(const G4UAdapter& rhs) const
{
  return (this == &rhs) ? true : false;
}

template <class UnplacedVolume_t>
G4UAdapter<UnplacedVolume_t>::
G4UAdapter(const G4UAdapter& rhs)
  : G4VSolid(rhs), UnplacedVolume_t(rhs),
    fRebuildPolyhedron(false), fPolyhedron(0)
{
  kHalfTolerance = 0.5*kCarTolerance;
}

template <class UnplacedVolume_t>
G4UAdapter<UnplacedVolume_t>& G4UAdapter<UnplacedVolume_t>::
operator=(const G4UAdapter& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fRebuildPolyhedron = false;
  delete fPolyhedron; fPolyhedron = 0;
  kHalfTolerance = 0.5*kCarTolerance;

  return *this;
}

template <class UnplacedVolume_t>
EInside G4UAdapter<UnplacedVolume_t>::
Inside(const G4ThreeVector& p) const
{
  U3Vector pt(p.x(), p.y(), p.z());
  vecgeom::EnumInside in_temp;
  EInside in = kOutside;

  in_temp = UnplacedVolume_t::Inside(pt);

  if (in_temp == vecgeom::EnumInside::eInside) in = kInside;
  else if (in_temp == vecgeom::EnumInside::eSurface) in = kSurface;

  return in;
}

template <class UnplacedVolume_t>
G4ThreeVector G4UAdapter<UnplacedVolume_t>::
SurfaceNormal(const G4ThreeVector& pt) const
{
  U3Vector p(pt.x(), pt.y(), pt.z());
  U3Vector n;
  UnplacedVolume_t::Normal(p, n);
  return G4ThreeVector(n.x(), n.y(), n.z());
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::
DistanceToIn(const G4ThreeVector& pt, const G4ThreeVector& d) const
{
  U3Vector p(pt.x(), pt.y(), pt.z());
  U3Vector v(d.x(), d.y(), d.z());
  G4double dist = UnplacedVolume_t::DistanceToIn(p, v, kInfinity);

  // apply Geant4 distance conventions
  //
  if (dist < kHalfTolerance)  return 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::
DistanceToIn(const G4ThreeVector& pt) const
{
  U3Vector p(pt.x(), pt.y(), pt.z());
  G4double dist = UnplacedVolume_t::SafetyToIn(p);

  // Apply Geant4 convention: convert negative values to zero
  //
  if (dist < kHalfTolerance)  return 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::
DistanceToOut(const G4ThreeVector& pt, const G4ThreeVector& d,
              const G4bool calcNorm, G4bool* validNorm,
                    G4ThreeVector* norm) const
{
  U3Vector p(pt.x(), pt.y(), pt.z());
  U3Vector v(d.x(), d.y(), d.z());

  G4double dist = UnplacedVolume_t::DistanceToOut(p, v, kInfinity);
  if(calcNorm) // *norm=n, but only after calcNorm check and if convex volume
  {
    if (UnplacedVolume_t::IsConvex())
    {
      U3Vector n, hitpoint = p + dist * v;
      UnplacedVolume_t::Normal(hitpoint, n);
      *validNorm = true;
      norm->set(n.x(), n.y(), n.z());
    }
    else
    {
      *validNorm = false;
    }
  }

  // Apply Geant4 distance conventions
  //
  if (dist < kHalfTolerance)  return 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::
DistanceToOut(const G4ThreeVector& pt) const
{
  U3Vector p(pt.x(), pt.y(), pt.z());
  G4double dist = UnplacedVolume_t::SafetyToOut(p);

  // Apply Geant4 convention: convert negative values to zero
  //
  if (dist < kHalfTolerance)  return 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::GetCubicVolume()
{
  return UnplacedVolume_t::Capacity();
}

template <class UnplacedVolume_t>
G4double G4UAdapter<UnplacedVolume_t>::GetSurfaceArea()
{
  return UnplacedVolume_t::SurfaceArea();
}

template <class UnplacedVolume_t>
G4ThreeVector G4UAdapter<UnplacedVolume_t>::GetPointOnSurface() const
{
  U3Vector p = UnplacedVolume_t::SamplePointOnSurface();;
  return G4ThreeVector(p.x(), p.y(), p.z());
}

// Inline visualization adapters

namespace
{
  G4Mutex pMutex = G4MUTEX_INITIALIZER;
}

template <class UnplacedVolume_t>
void G4UAdapter<UnplacedVolume_t>::
ComputeDimensions(G4VPVParameterisation*, const G4int,
                                          const G4VPhysicalVolume*)
{
    std::ostringstream message;
    message << "Illegal call to G4UAdapter::ComputeDimensions()" << G4endl
            << "Method not overloaded by derived class !";
    G4Exception("G4UAdapter::ComputeDimensions()", "GeomSolids0003",
                FatalException, message);
}

template <class UnplacedVolume_t>
void G4UAdapter<UnplacedVolume_t>::
DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

template <class UnplacedVolume_t>
G4GeometryType G4UAdapter<UnplacedVolume_t>::
GetEntityType() const
{

  G4String string = "VSolid"; // UnplacedVolume_t::GetEntityType();
  return "G4" + string;
}

template <class UnplacedVolume_t>
std::ostream& G4UAdapter<UnplacedVolume_t>::
StreamInfo(std::ostream& os) const
{
  UnplacedVolume_t::Print(os);
  return os;
}

template <class UnplacedVolume_t>
G4VSolid* G4UAdapter<UnplacedVolume_t>::Clone() const
{
  std::ostringstream message;
  message << "Clone() method not implemented for type: "
          << GetEntityType() << "!" << G4endl
          << "Returning NULL pointer!";
  G4Exception("G4UAdapter::Clone()", "GeomSolids1001", JustWarning, message);
  return 0;
}

template <class UnplacedVolume_t>
G4bool G4UAdapter<UnplacedVolume_t>::CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const
{
  U3Vector vmin, vmax;
  UnplacedVolume_t::Extent(vmin,vmax);
  G4ThreeVector bmin(vmin.x(),vmin.y(),vmin.z());
  G4ThreeVector bmax(vmax.x(),vmax.y(),vmax.z());

  // Check correctness of the bounding box
  //
  if (bmin.x() >= bmax.x() || bmin.y() >= bmax.y() || bmin.z() >= bmax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " - " << GetEntityType() << " !"
            << "\nmin = " << bmin
            << "\nmax = " << bmax;
    G4Exception("G4UAdapter::CalculateExtent()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

template <class UnplacedVolume_t>
G4Polyhedron* G4UAdapter<UnplacedVolume_t>::CreatePolyhedron() const
{
  // Must be implemented in concrete wrappers...

  std::ostringstream message;
  message << "Visualization not supported for USolid shape "
          << GetEntityType() << "... Sorry!" << G4endl;
  G4Exception("G4UAdapter::CreatePolyhedron()", "GeomSolids0003",
              FatalException, message);
  return 0;
}

template <class UnplacedVolume_t>
G4Polyhedron* G4UAdapter<UnplacedVolume_t>::GetPolyhedron() const
{
  if (!fPolyhedron ||
      fRebuildPolyhedron ||
      fPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&pMutex);
    delete fPolyhedron;
    fPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fPolyhedron;
}

template <class UnplacedVolume_t>
G4VisExtent G4UAdapter<UnplacedVolume_t>::GetExtent() const
{
  U3Vector vmin, vmax;
  UnplacedVolume_t::Extent(vmin,vmax);
  return G4VisExtent(vmin.x(),vmax.x(),
                     vmin.y(),vmax.y(),
                     vmin.z(),vmax.z());
}

#endif  // G4GEOM_USE_USOLIDS

#endif // G4UADAPTER_HH
