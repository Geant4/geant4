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
// G4UAdapter
//
// Class description:
//
// Utility class for adapting VecGeom solids API to Geant4 solids.
// NOTE: Using protected inheritance since the Adapter is supposed to
// be a G4VSolid "implemented-in-terms-of" the VecGeom UnplacedVolume_t.
// The choice of protected vs private is due to the fact that we want
// to propagate functions further down in the inheritance hierarchy.

// Author: Gabriele Cosmo (CERN), 17.05.2017
//         Adapted for G4VSolid from original G4USolids bridge
//         class and the USolidsAdapter class in VecGeom.
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

#include "G4GeomTypes.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>

class G4VPVParameterisation;

/**
 * @brief G4UAdapter is a utility class for adapting VecGeom solids API to
 * Geant4 solids. The Adapter is supposed to be a G4VSolid
 * "implemented-in-terms-of" the VecGeom UnplacedVolume_t.
 */

template <class UnplacedVolume_t>
class G4UAdapter : public G4VSolid, protected UnplacedVolume_t
{
  public:

    using U3Vector = vecgeom::Vector3D<G4double>;

    /** VecGeom volumes have special delete/new ("AlignedBase")
        and we need to make these functions public again. */
    using UnplacedVolume_t::operator delete;
    using UnplacedVolume_t::operator new;

    /**
     * Constructor taking a name.
     *  @param[in] name The name of the volume.
     */
    G4UAdapter(const G4String& name);

    /**
     * Constructor templated on arguments for UnplacedVolume_t.
     *  @param[in] name The name of the volume.
     *  @param[in] params Templated arguments for UnplacedVolume_t.
     */
    template <typename... T>
    G4UAdapter(const G4String& name, const T &... params);

    /**
     * Virtual destructor.
     */
    virtual ~G4UAdapter();

    /**
     * Copy constructor and assignment operator.
     */
    G4UAdapter(const G4UAdapter& rhs);
    G4UAdapter& operator=(const G4UAdapter& rhs);

    /**
     * Equality operator. Returns true only if addresses are the same.
     */
    G4bool operator==(const G4UAdapter& s) const;

    /**
     * Calculates the minimum and maximum extent of the solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                               G4double& pMin, G4double& pMax) const override;

    /**
     * Returns the characterisation of a point at offset 'p' respect
     * to the shape.
     *  @param[in] p The point at offset p.
     *  @returns kOutside if the point is outside the shapes boundaries
     *           plus Tolerance/2; kSurface if the point is less than
     *           Tolerance/2 from a surface; kInside otherwise.
     */
    virtual EInside Inside(const G4ThreeVector& p) const override;

    /**
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset 'p'.
     *  @param[in] p The point at offset p.
     *  @returns The outwards pointing unit normal.
     */
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const override;

    /**
     * Returns the distance along the normalised vector 'v' to the shape,
     * from the point at offset 'p'. If there is no intersection, returns
     * kInfinity. The first intersection resulting from 'leaving' a
     * surface/volume is discarded. Hence, it is tolerant of points on
     * the surface of the shape.
     *  @param[in] p The point at offset p.
     *  @param[in] v The normalised direction vector.
     *  @returns The distance to enter the shape.
     */
    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v) const override;

    /**
     * Calculates the distance to the nearest surface of a shape from an
     * outside point. The distance can be an underestimate.
     *  @param[in] p The point at offset p.
     *  @returns The safety distance to enter the shape.
     */
    virtual G4double DistanceToIn(const G4ThreeVector& p) const override;

    /**
     * Returns the distance along the normalised vector 'v' to the shape,
     * from a point at an offset 'p' inside or on the surface of the shape.
     * Intersections with surfaces, when the point is less than Tolerance/2
     * from a surface must be ignored.
     *  @param[in] p The point at offset p.
     *  @param[in] v The normalised direction vector.
     *  @param[in] calcNorm Flag to indicate if to calculate the normal or not.
     *  @param[out] validNorm Flag set to true if the solid lies entirely
     *              behind or on the exiting surface. It is set false if the
     *              solid does not lie entirely behind or on the exiting surface.
     *              'calcNorm' must be true, otherwise it is unused.
     *  @param[out] n The exiting outwards normal vector (undefined Magnitude).
     *              'calcNorm' must be true, otherwise it is unused.
     *  @returns The distance to exit the shape.
     */
    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm = false,
                                   G4bool* validNorm = 0,
                                   G4ThreeVector* n = 0) const override;

    /**
     * Calculates the distance to the nearest surface of a shape from an
     * inside point 'p'. The distance can be an underestimate.
     *  @param[in] p The point at offset p.
     *  @returns The safety distance to exit the shape.
     */
    virtual G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation. Throws exception if ComputeDimensions() is
     * called from an illegal derived class.
     */
    virtual void ComputeDimensions(G4VPVParameterisation* p,
                                   const G4int n,
                                   const G4VPhysicalVolume* pRep) override;

    /**
     * Returns an estimation of the solid volume in internal units.
     * This method may be overloaded by derived classes to compute the
     * exact geometrical quantity for solids where this is possible,
     * or anyway to cache the computed value.
     * Note: the computed value is NOT cached.
     */
    virtual G4double GetCubicVolume() override;

    /**
     * Returns an estimation of the solid surface area in internal units.
     * This method may be overloaded by derived classes to compute the
     * exact geometrical quantity for solids where this is possible,
     * or anyway to cache the computed value.
     * Note: the computed value is NOT cached.
     */
    virtual G4double GetSurfaceArea() override;

    /**
     * Returns a random point located on the surface of the solid.
     */
    virtual G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns the number of constituents used for construction of the solid.
     * For non-Boolean solids the return value is one.
     */
    virtual G4int GetNumOfConstituents() const override;

    /**
     * Returns true if the solid has only planar faces, false otherwise.
     */
    virtual G4bool IsFaceted() const override;

    /**
     * Provides identification of the class of an object
     * (required for persistency).
     */
    virtual G4GeometryType GetEntityType() const override;

    /**
     * Returns a pointer of a dynamically allocated copy of the solid.
     * Returns a null pointer with warning in case the concrete solid does not
     * implement this method. The caller has responsibility for ownership.
     */
    virtual G4VSolid* Clone() const override;

    /**
     * Dumps contents of the solid to a stream.
     */
    virtual std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * A "double dispatch" function which identifies the solid
     * to the graphics scene for visualization.
     */
    virtual void DescribeYourselfTo(G4VGraphicsScene& scene) const override;

    /**
     * Provides extent (bounding box) as possible hint to the graphics view.
     */
    virtual G4VisExtent GetExtent()  const override;

    /**
     * Creates a Polyhedron used for Visualisation.
     */
    virtual G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Smart access function - creates on request and stores for future
     * access. A null pointer means "not available".
     */
    virtual G4Polyhedron* GetPolyhedron() const override;


    // VecGeom overridden methods ---------------------------------------------

    vecgeom::Precision
    DistanceToOut(U3Vector const& position, U3Vector const& direction,
                  vecgeom::Precision stepMax = kInfinity) const override
    {
      return UnplacedVolume_t::DistanceToOut(position, direction, stepMax);
    }

    vecgeom::EnumInside
    Inside(U3Vector const& aPoint) const override
    {
      return UnplacedVolume_t::Inside(aPoint);
    }

    vecgeom::Precision
    DistanceToIn(U3Vector const& position, U3Vector const& direction,
                 const vecgeom::Precision step_max = kInfinity) const override
    {
      return UnplacedVolume_t::DistanceToIn(position, direction, step_max);
    }

    G4bool Normal(U3Vector const& aPoint, U3Vector& aNormal) const override
    {
      return UnplacedVolume_t::Normal(aPoint, aNormal);
    }

    void Extent(U3Vector& aMin, U3Vector& aMax) const override
    {
      return UnplacedVolume_t::Extent(aMin, aMax);
    }

    U3Vector SamplePointOnSurface() const override
    {
      return UnplacedVolume_t::SamplePointOnSurface();
    }

  protected:  // data

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fPolyhedron = nullptr;

    /** Cached geometrical tolerance. */
    G4double kHalfTolerance;

    using UnplacedVolume_t::DistanceToOut;
    using UnplacedVolume_t::DistanceToIn;
};

// Inline implementations

#include "G4UAdapter.icc"

#endif  // G4GEOM_USE_USOLIDS

#endif // G4UADAPTER_HH
