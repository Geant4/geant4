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
// G4Polycone
//
// Class description:
//
// Class implementing a CSG-like type "PCON" Geant 3.21 volume,
// inherited from  class G4VCSGfaceted:
//
//   G4Polycone( const G4String& name,
//               G4double phiStart,     // initial phi starting angle
//               G4double phiTotal,     // total phi angle
//               G4int numZPlanes,      // number of z planes
//               const G4double zPlane[],  // position of z planes
//               const G4double rInner[],  // tangent distance to inner surface
//               const G4double rOuter[])  // tangent distance to outer surface
//
// Alternative constructor:
//
//   G4Polycone( const G4String& name,
//               G4double phiStart,   // initial phi starting angle
//               G4double phiTotal,   // total phi angle
//               G4int numRZ,         // number corners in r,z space
//               const G4double r[],  // r coordinates of these corners
//               const G4double z[])  // z coordinates of these corners

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4POLYCONE_HH
#define G4POLYCONE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPOLYCONE 1
#endif

#if defined(G4GEOM_USE_UPOLYCONE)
  #define G4UPolycone G4Polycone
  #include "G4UPolycone.hh"
#else

#include "G4VCSGfaceted.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyconeHistorical.hh"
#include "G4Polyhedron.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4VCSGface;

/**
 * @brief G4Polycone represents a composed closed shape (PCON) made of
 * cones and cylinders, along the Z axis with increasing Z, with or without
 * cut in Phi.
 */

class G4Polycone : public G4VCSGfaceted
{
  public:

    /**
     * Constructs a polycone shape, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numZPlanes Number of Z planes.
     *  @param[in] zPlane Position of Z planes, with Z in increasing order.
     *  @param[in] rInner Tangent distance to inner surface.
     *  @param[in] rOuter Tangent distance to outer surface.
     */
    G4Polycone( const G4String& name,
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int numZPlanes,     // number of z planes
                const G4double zPlane[],    // position of z planes
                const G4double rInner[],    // tangent distance to inner surface
                const G4double rOuter[]  ); // tangent distance to outer surface

    /**
     * Alternative constructor of a polycone shape, given corners coordinates.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numRZ Number of corners in r,Z space.
     *  @param[in] r r coordinates of corners.
     *  @param[in] z Z coordinates of corners.
     */
    G4Polycone( const G4String& name,
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int numRZ,          // number corners in r,z space
                const G4double r[],         // r coordinates of these corners
                const G4double z[]       ); // z coordinates of these corners

    /**
     * Destructor.
     */
    ~G4Polycone() override;

    /**
     * Accessors.
     */
    inline G4double GetStartPhi()    const;
    inline G4double GetEndPhi()      const;
    inline G4double GetSinStartPhi() const;
    inline G4double GetCosStartPhi() const;
    inline G4double GetSinEndPhi()   const;
    inline G4double GetCosEndPhi()   const;
    inline G4bool IsOpen()           const;
    inline G4int GetNumRZCorner()    const;
    inline G4PolyconeSideRZ GetCorner(G4int index) const;

    /**
     * Gets and sets the original parameters of the solid.
     */
    inline G4PolyconeHistorical* GetOriginalParameters() const;
    inline void SetOriginalParameters(G4PolyconeHistorical* pars);

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in G4VSolid. Remaining functions are concretely
     * defined in the base class G4VCSGfaceted.
     */
    EInside Inside( const G4ThreeVector& p ) const override;
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn( const G4ThreeVector& p ) const override;

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
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
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions( G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    /**
     * Returns the type ID, "G4Polycone" of the solid.
     */
    G4GeometryType GetEntityType() const override;

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
     * Returns a pointer to a generated polyhedron used for visualisation.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Clears all parameters and rebuild the shape, for use in divisions.
     */
    G4bool Reset();

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Polycone(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Polycone( const G4Polycone& source );
    G4Polycone& operator=( const G4Polycone& source );

  private:

    /**
     * Generic initializer, called by all constructors.
     */
    G4bool SetOriginalParameters(G4ReduciblePolygon* rz);

    void Create( G4double phiStart,        // initial phi starting angle
                 G4double phiTotal,        // total phi angle
                 G4ReduciblePolygon* rz ); // r/z coordinate of these corners

    /**
     * Copy parameters from other solid; used in copy constructor and
     * assignment operator.
     */
    void CopyStuff( const G4Polycone& source );

    /**
     * Sets the vector of surface elements. Auxiliary method used for
     * sampling random points on surface.
     */
    void SetSurfaceElements() const;

  private:

    /** The original parameters. */
    G4double startPhi;        // Starting phi value (0 < phiStart < 2pi)
    G4double endPhi;          // End phi value (0 < endPhi-phiStart < 2pi)
    G4bool phiIsOpen = false; // True if there is a phi segment
    G4int numCorner;          // Number RZ points
    G4PolyconeSideRZ* corners = nullptr; // Corner r,z points
    G4PolyconeHistorical* original_parameters = nullptr; // Original input pars

    G4EnclosingCylinder* enclosingCylinder = nullptr; // Our quick test

    struct surface_element { G4double area = 0.; G4int i0 = 0, i1 = 0, i2 = 0; };
    mutable std::vector<surface_element>* fElements = nullptr;
};

#include "G4Polycone.icc"

#endif

#endif
