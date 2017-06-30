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
// $Id:$
// GEANT4 tag $Name:$
//
//
// class G4USolid
//
// Class description:
//
// Bridge base class for solids defined in the Unified Solids Library.

// --------------------------------------------------------------------
#ifndef G4USolid_HH
#define G4USolid_HH

#include "G4VSolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "VUSolid.hh"

class G4VPVParameterisation;

class G4USolid : public G4VSolid
{
  public:  // with description

    G4USolid(const G4String& pName, VUSolid* shape);
    // Creates a new shape, with the supplied name. No provision is made
    // for sharing a common name amongst multiple classes.
    virtual ~G4USolid();
    // Default destructor.

    G4bool operator==(const G4USolid& s) const;
    // Return true only if addresses are the same.

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                   G4double& pMin, G4double& pMax) const;
    // Calculate the minimum and maximum extent of the solid, when under the
    // specified transform, and within the specified limits. If the solid
    // is not intersected by the region, return false, else return true.

    virtual EInside Inside(const G4ThreeVector& p) const;
    // Returns kOutside if the point at offset p is outside the shapes
    // boundaries plus Tolerance/2, kSurface if the point is <= Tolerance/2
    // from a surface, otherwise kInside.

    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
    // Returns the outwards pointing unit normal of the shape for the
    // surface closest to the point at offset p.

    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const;
    // Return the distance along the normalised vector v to the shape,
    // from the point at offset p. If there is no intersection, return
    // kInfinity. The first intersection resulting from `leaving' a
    // surface/volume is discarded. Hence, it is tolerant of points on
    // the surface of the shape.

    virtual G4double DistanceToIn(const G4ThreeVector& p) const;
    // Calculate the distance to the nearest surface of a shape from an
    // outside point. The distance can be an underestimate.

    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm = false,
                                   G4bool* validNorm = 0,
                                   G4ThreeVector* n = 0) const;
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

    virtual G4double DistanceToOut(const G4ThreeVector& p) const;
    // Calculate the distance to the nearest surface of a shape from an
    // inside point. The distance can be an underestimate.

    virtual void ComputeDimensions(G4VPVParameterisation* p,
                                   const G4int n,
                                   const G4VPhysicalVolume* pRep);
      // Throw exception if ComputeDimensions called from an illegal
      // derived class.

    virtual G4double GetCubicVolume();
    // Returns an estimation of the solid volume in internal units.
    // This method may be overloaded by derived classes to compute the
    // exact geometrical quantity for solids where this is possible,
    // or anyway to cache the computed value.
    // Note: the computed value is NOT cached.

    virtual G4double GetSurfaceArea();
    // Return an estimation of the solid surface area in internal units.
    // This method may be overloaded by derived classes to compute the
    // exact geometrical quantity for solids where this is possible,
    // or anyway to cache the computed value.
    // Note: the computed value is NOT cached.

    virtual G4GeometryType  GetEntityType() const;
    // Provide identification of the class of an object.
    // (required for persistency and STEP interface)

    virtual G4ThreeVector GetPointOnSurface() const;
    // Returns a random point located on the surface of the solid.

    virtual G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.
    // Returns NULL pointer with warning in case the concrete solid does not
    // implement this method. The caller has responsibility for ownership.

    virtual std::ostream& StreamInfo(std::ostream& os) const;
    // Dumps contents of the solid to a stream.

    virtual void DescribeYourselfTo(G4VGraphicsScene& scene) const;
    // A "double dispatch" function which identifies the solid
    // to the graphics scene for visualization.

    virtual G4VisExtent   GetExtent()  const;
    // Provide extent (bounding box) as possible hint to the graphics view.
    virtual G4Polyhedron* CreatePolyhedron() const;
    // Create Polyhedron used for Visualisation
    virtual G4Polyhedron* GetPolyhedron() const;
    // Smart access function - creates on request and stores for future
    // access.  A null pointer means "not available".

  public:  // without description

    G4USolid(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    G4USolid(const G4USolid& rhs);
    G4USolid& operator=(const G4USolid& rhs);
    // Copy constructor and assignment operator.

    VUSolid* GetSolid() const
    {
      return fShape;
    }

  protected:  // data

    VUSolid* fShape;
    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fPolyhedron;
};

// Inline implementations

inline
EInside G4USolid::Inside(const G4ThreeVector& p) const
{
  UVector3 pt;
  VUSolid::EnumInside in_temp;
  EInside in = kOutside;
  pt.x() = p.x();
  pt.y() = p.y();
  pt.z() = p.z(); // better assign at construction

  in_temp = fShape->Inside(pt);

  if (in_temp == VUSolid::EnumInside::eInside) in = kInside;
  else if (in_temp == VUSolid::EnumInside::eSurface) in = kSurface;

  return in;
}

inline
G4ThreeVector G4USolid::SurfaceNormal(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z();
  UVector3 n;
  fShape->Normal(p, n);
  return G4ThreeVector(n.x(), n.y(), n.z());
}

inline
G4double G4USolid::DistanceToIn(const G4ThreeVector& pt,
                                const G4ThreeVector& d) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  UVector3 v;
  v.x() = d.x();
  v.y() = d.y();
  v.z() = d.z(); // better assign at construction
  G4double dist = fShape->DistanceToIn(p, v);
//  return  (dist > halfTolerance) ? dist : 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

inline
G4double G4USolid::DistanceToIn(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  G4double dist = fShape->SafetyFromOutside(p); // true?
//  return (dist > halfTolerance) ? dist : 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

inline
G4double G4USolid::DistanceToOut(const G4ThreeVector& pt,
                                 const G4ThreeVector& d,
                                 const G4bool calcNorm,
                                 G4bool* validNorm,
                                 G4ThreeVector* norm) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  UVector3 v;
  v.x() = d.x();
  v.y() = d.y();
  v.z() = d.z(); // better assign at construction
  UVector3 n;
  G4bool valid;
  G4double dist = fShape->DistanceToOut(p, v, n, valid); // should use local variable
  if(calcNorm)
  {
    if(valid){ *validNorm = true; }
    else { *validNorm = false; }
    if(*validNorm)  // *norm = n, but only after calcNorm check
    {
      norm->setX(n.x());
      norm->setY(n.y());
      norm->setZ(n.z());
    }
  }
//  return (dist > halfTolerance) ? dist : 0.0;
  return (dist > kInfinity) ? kInfinity : dist;
}

inline
G4double G4USolid::DistanceToOut(const G4ThreeVector& pt) const
{
  UVector3 p;
  p.x() = pt.x();
  p.y() = pt.y();
  p.z() = pt.z(); // better assign at construction
  G4double dist = fShape->SafetyFromInside(p); // true?
//  return (dist > halfTolerance) ? dist : 0.0;
  return dist;
}

inline
G4double G4USolid::GetCubicVolume()
{
  return fShape->Capacity();
}

inline
G4double G4USolid::GetSurfaceArea()
{
  return fShape->SurfaceArea();
}

inline
G4ThreeVector G4USolid::GetPointOnSurface() const
{
  UVector3 p;
  p = fShape->GetPointOnSurface();
  return G4ThreeVector(p.x(), p.y(), p.z());
}

#endif  // G4GEOM_USE_USOLIDS

#endif
