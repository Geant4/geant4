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
// G4GenericPolycone
//
// Class description:
//
// Implementing a GenericPolycone constructed by points with (r,z)
// coordinates and allowing Z 'go back'.
//
//   G4GenericPolycone( const G4String& name,
//                      G4double phiStart,   // initial phi starting angle
//                      G4double phiTotal,   // total phi angle
//                      G4int    numRZ,      // number corners in r,z space
//                      const G4double r[],  // r coordinate of these corners
//                      const G4double z[])  // z coordinate of these corners

// Authors: T.Nikitina, G.Cosmo (CERN), 29.10.2013 - Created
// --------------------------------------------------------------------
#ifndef G4GENERICPOLYCONE_HH
#define G4GENERICPOLYCONE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UGENERICPOLYCONE 1
#endif

#if defined(G4GEOM_USE_UGENERICPOLYCONE)
  #define G4UGenericPolycone G4GenericPolycone
  #include "G4UGenericPolycone.hh"
#else

#include "G4VCSGfaceted.hh"
#include "G4PolyconeSide.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4VCSGface;

/**
 * @brief G4GenericPolycone is a Polycone shape where the composing Z planes
 * positions, in their order of definition, may not be monotically increasing,
 * i.e. may also decrease.
 */

class G4GenericPolycone : public G4VCSGfaceted
{
  public:

    /**
     * Constructs a generic polycone shape, given its parameters.
     *  @param[in] name The solid name.
     *  @param[in] phiStart The initial Phi starting angle.
     *  @param[in] phiTotal The total Phi angle.
     *  @param[in] numRZ Number of corners in r,Z space.
     *  @param[in] r Vector of r coordinate of corners.
     *  @param[in] z Vector of Z coordinate of corners.
     */
    G4GenericPolycone( const G4String& name,
                             G4double phiStart,
                             G4double phiTotal,
                             G4int    numRZ,
                       const G4double r[],
                       const G4double z[] );

    /**
     * Destructor.
     */
    ~G4GenericPolycone() override;

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
    inline G4int  GetNumRZCorner()   const;
    inline G4PolyconeSideRZ GetCorner(G4int index) const;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside( const G4ThreeVector &p ) const override;
    G4double DistanceToIn( const G4ThreeVector &p,
                           const G4ThreeVector &v ) const override;
    G4double DistanceToIn( const G4ThreeVector &p ) const override;
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;
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
     * Returns the type ID, "G4GenericPolycone" of the solid.
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
     * Returns a pointer to a polyhedron for use in visualisation.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Does nothing. Reset of parameters (for use in divisions) is not
     * allowed for a generic polycone. Issues a warning and just returns true.
     */
    G4bool Reset();

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4GenericPolycone(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4GenericPolycone( const G4GenericPolycone& source );
    G4GenericPolycone& operator=( const G4GenericPolycone& source );

  private:

    /**
     * Generic initializer, called by constructor.
     */
    void Create( G4double phiStart,        // initial phi starting angle
                 G4double phiTotal,        // total phi angle
                 G4ReduciblePolygon* rz ); // r/z coordinate of these corners

    /**
     * Utility for copying contents, used in copy constructor and assignment
     * operator.
     */
    void CopyStuff( const G4GenericPolycone& source );

     /**
      * Auxiliary method for sampling random points on surface.
      * Sets the vector of surface elements.
      */
    void SetSurfaceElements() const;

  private:

    /** Original parameters. */
    G4double startPhi;            // Starting phi value (0 < phiStart < 2pi)
    G4double endPhi;              // end phi value (0 < endPhi-phiStart < 2pi)
    G4bool   phiIsOpen = false;   // true if there is a phi segment
    G4int    numCorner;           // number RZ points
    G4PolyconeSideRZ* corners = nullptr;  // corner r,z points

    G4EnclosingCylinder* enclosingCylinder = nullptr; // Our quick test

    struct surface_element { G4double area = 0.; G4int i0 = 0, i1 = 0, i2 = 0; };
    mutable std::vector<surface_element>* fElements = nullptr;
};

#include "G4GenericPolycone.icc"

#endif

#endif
