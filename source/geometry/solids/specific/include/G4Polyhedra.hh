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
// G4Polyhedra
//
// Class description:
//
// Class implementing a CSG-like type "PGON" Geant 3.21 volume,
// inherited from class G4VCSGfaceted:
//
//   G4Polyhedra( const G4String& name,
//                G4double phiStart,         - initial phi starting angle
//                G4double phiTotal,         - total phi angle
//                G4int numSide,             - number sides
//                G4int numZPlanes,          - number of z planes
//                const G4double zPlane[],   - position of z planes
//                const G4double rInner[],   - tangent distance to inner surface
//                const G4double rOuter[]  ) - tangent distance to outer surface
//
//   G4Polyhedra( const G4String& name,
//                G4double phiStart,    - initial phi starting angle
//                G4double phiTotal,    - total phi angle
//                G4int    numSide,     - number sides
//                G4int    numRZ,       - number corners in r,z space
//                const G4double r[],   - r coordinates of these corners
//                const G4double z[] )  - z coordinates of these corners

// Author: David C. Williams (UCSC), 1998 - First implementation
// --------------------------------------------------------------------
#ifndef G4POLYHEDRA_HH
#define G4POLYHEDRA_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPOLYHEDRA 1
#endif

#if defined(G4GEOM_USE_UPOLYHEDRA)
  #define G4UPolyhedra G4Polyhedra
  #include "G4UPolyhedra.hh"
#else

#include "G4VCSGfaceted.hh"
#include "G4PolyhedraSide.hh"
#include "G4PolyhedraHistorical.hh"
#include "G4Polyhedron.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;

/**
 * @brief G4Polyhedra represents a composed closed polyhedra (PGON) made of
 * planar sizes along the Z axis, with or without cut in Phi.
 */

class G4Polyhedra : public G4VCSGfaceted
{
  public:

    /**
     * Constructs a polyhedra, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numSide Number of sides.
     *  @param[in] numZPlanes Number of Z planes.
     *  @param[in] zPlane Position of Z planes.
     *  @param[in] rInner Tangent distance to inner surface.
     *  @param[in] rOuter Tangent distance to outer surface.
     */
    G4Polyhedra(const G4String& name,
                      G4double phiStart,   // initial phi starting angle
                      G4double phiTotal,   // total phi angle
                      G4int numSide,       // number sides
                      G4int numZPlanes,    // number of z planes
                const G4double zPlane[],   // position of z planes
                const G4double rInner[],   // tangent distance to inner surface
                const G4double rOuter[] ); // tangent distance to outer surface

    /**
     * Alternative constructor of a polyhedra, given corners coordinates.
     *  @param[in] name The solid name.
     *  @param[in] phiStart Initial Phi starting angle.
     *  @param[in] phiTotal Total Phi angle.
     *  @param[in] numSide Number of sides.
     *  @param[in] numRZ Number of corners in r,Z space.
     *  @param[in] r r coordinates of corners.
     *  @param[in] z Z coordinates of corners.
     */
    G4Polyhedra(const G4String& name,
                      G4double phiStart,   // initial phi starting angle
                      G4double phiTotal,   // total phi angle
                      G4int    numSide,    // number sides
                      G4int    numRZ,      // number corners in r,z space
                const G4double r[],        // r coordinates of these corners
                const G4double z[] );      // z coordinates of these corners

    /**
     * Destructor.
     */
    ~G4Polyhedra() override;

    /**
     * Accessors.
     */
    inline G4int GetNumSide()        const;
    inline G4double GetStartPhi()    const;
    inline G4double GetEndPhi()      const;
    inline G4double GetSinStartPhi() const;
    inline G4double GetCosStartPhi() const;
    inline G4double GetSinEndPhi()   const;
    inline G4double GetCosEndPhi()   const;
    inline G4bool IsOpen()           const;
    inline G4bool IsGeneric()        const;
    inline G4int GetNumRZCorner()    const;
    inline G4PolyhedraSideRZ GetCorner( const G4int index ) const;

    /**
     * Returns internal scaled parameters.
     */
    inline G4PolyhedraHistorical* GetOriginalParameters() const;

    /**
     * Sets internal parameters. Parameters 'Rmin' and 'Rmax' in input must
     * be scaled first by a factor computed as 'cos(0.5*phiTotal/theNumSide)',
     * if not already scaled.
     */
    inline void SetOriginalParameters(G4PolyhedraHistorical* pars);

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
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    /**
     * Returns the type ID, "G4Polyhedra" of the solid.
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
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo( std::ostream& os ) const override;

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
    G4Polyhedra(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Polyhedra( const G4Polyhedra& source );
    G4Polyhedra& operator=( const G4Polyhedra& source );

  private:

    /**
     * Sets internal parameters for the generic constructor.
     */
    void SetOriginalParameters(G4ReduciblePolygon* rz);

    /**
     * Generates the shape and is called by each constructor,
     * after the conversion of the arguments.
     */
    void Create( G4double phiStart,           // initial phi starting angle
                 G4double phiTotal,           // total phi angle
                 G4int    numSide,            // number sides
                 G4ReduciblePolygon* rz );    // rz coordinates

    /**
     * Copy parameters from other solid or reset them.
     * Used in copy constructor and assignment operator.
     */
    void CopyStuff( const G4Polyhedra& source );
    void DeleteStuff();

    /**
     * Sets the vector of surface elements. Auxiliary method used for
     * sampling random points on surface.
     */
    void SetSurfaceElements() const;

  private:

    G4int numSide = 0;    // Number of sides
    G4double startPhi;    // Starting phi value (0 < phiStart < 2pi)
    G4double endPhi;      // end phi value (0 < endPhi-phiStart < 2pi)
    G4bool phiIsOpen = false;   // true if there is a phi segment
    G4bool genericPgon = false; // true if created through 2nd generic ctor
    G4int numCorner = 0;  // number RZ points
    G4PolyhedraSideRZ* corners = nullptr;  // our corners
    G4PolyhedraHistorical* original_parameters = nullptr; // original input pars

    G4EnclosingCylinder* enclosingCylinder = nullptr;

    struct surface_element { G4double area=0.; G4int i0=0, i1=0, i2=0; };
    mutable std::vector<surface_element>* fElements = nullptr;
};

#include "G4Polyhedra.icc"

#endif

#endif
