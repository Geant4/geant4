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
// * and NASA under contract number NNG04CT05P.                       *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// G4Tet
//
// Class description:
//
// A G4Tet is a tetrahedra solid, defined by 4 points in space.

// Author: M.H.Mendenhall & R.A.Weller (Vanderbilt University, USA), 03.09.2004
//         E.Tcherniaev (CERN), 08.01.2020 - Complete revision, speed up
// --------------------------------------------------------------------
#ifndef G4TET_HH
#define G4TET_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTET 1
#endif

#if defined(G4GEOM_USE_UTET)
  #define G4UTet G4Tet
  #include "G4UTet.hh"
#else

#include "G4VSolid.hh"

/**
 * @brief G4Tet is a tetrahedra solid, defined by 4 points in space.
 */

class G4Tet : public G4VSolid
{
  public:

    /**
     * Constructs a tetrahedra, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] anchor The anchor point.
     *  @param[in] p2 Point 2.
     *  @param[in] p3 Point 3.
     *  @param[in] p4 Point 4.
     *  @param[in] degeneracyFlag Flag indicating degeneracy of points.
     */
    G4Tet(const G4String& pName,
          const G4ThreeVector& anchor,
          const G4ThreeVector& p2,
          const G4ThreeVector& p3,
          const G4ThreeVector& p4,
                G4bool* degeneracyFlag = nullptr);

    /**
     * Destructor.
     */
    ~G4Tet() override;

    /**
     * Modifier and accessors, for the four vertices of the shape.
     */
    void SetVertices(const G4ThreeVector& anchor,
                     const G4ThreeVector& p1,
                     const G4ThreeVector& p2,
                     const G4ThreeVector& p3,
                     G4bool* degeneracyFlag = nullptr);
    void GetVertices(G4ThreeVector& anchor,
                     G4ThreeVector& p1,
                     G4ThreeVector& p2,
                     G4ThreeVector& p3) const;
    std::vector<G4ThreeVector> GetVertices() const;

    /**
     * Checks if the tetrahedron is degenerate. A tetrahedron is considered
     * as degenerate in case its minimal height is less than the degeneracy
     * tolerance
     *  @returns true if the tetrahedron is degenerate.
     */
    G4bool CheckDegeneracy(const G4ThreeVector& p0,
                           const G4ThreeVector& p1,
                           const G4ThreeVector& p2,
                           const G4ThreeVector& p3) const;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    void SetBoundingLimits(const G4ThreeVector& pMin, const G4ThreeVector& pMax);
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

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
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const override;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Tet" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    G4bool IsFaceted () const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent () const override;
    G4Polyhedron* CreatePolyhedron () const override;
    G4Polyhedron* GetPolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Tet(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Tet(const G4Tet& rhs);
    G4Tet& operator=(const G4Tet& rhs);

  private:

    /**
     * Initialises the data members.
     */
    void Initialize(const G4ThreeVector& p0,
                    const G4ThreeVector& p1,
                    const G4ThreeVector& p2,
                    const G4ThreeVector& p3);

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

  private:

    G4double halfTolerance = 0;
    G4double fCubicVolume = 0; // Volume
    G4double fSurfaceArea = 0; // Surface area
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

    G4ThreeVector fVertex[4];   // thetrahedron vertices
    G4ThreeVector fNormal[4];   // normals to faces
    G4double fDist[4] = {0};    // distances from origin to faces
    G4double fArea[4] = {0};    // face areas
    G4ThreeVector fBmin, fBmax; // bounding box
};

#endif

#endif
