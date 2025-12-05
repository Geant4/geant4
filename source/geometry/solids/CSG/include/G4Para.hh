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
// G4Para
//
// Class description:
//
// A parallelepiped, essentially a box with half lengths dx,dy,dz
// 'skewed' so that there are angles theta & phi of the polar line
// joining the faces at +-dz in z, and alpha formed by the y axis
// and the plane joining the centre of the faces parallel to the
// z-x plane at -dy and +dy.
//
// A G4Para is defined by:
//   dx,dy,dz - Half-length in x,y,z
//   alpha    - Angle formed by the y axis and by the plane joining
//              the centre of the faces parallel to the z-x plane
//              at -dy and +dy
//   theta    - Polar angle of the line joining the centres of the
//              faces at -dz and +dz in z
//   phi      - Azimuthal angle of the line joining the centres of the
//              faces at -dz and +dz in z
// Member data:
//
//   Note that the angles parameters are not stored - precomputed trig is
//   stored instead.
//
//      fDx   Half-length in x
//      fDy   Half-length in y
//      fDz   Half-length in z
//
//      fTalpha       Tan of alpha
//      fTthetaCphi   Tan theta * Cos phi
//      fTthetaSphi   Tan theta * Sin phi

// Author: Paul Kent (CERN), 21.03.1994 - Code converted to tolerant geometry
// --------------------------------------------------------------------
#ifndef G4PARA_HH
#define G4PARA_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPARA 1
#endif

#if defined(G4GEOM_USE_UPARA)
  #define G4UPara G4Para
  #include "G4UPara.hh"
#else

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

/**
 * @brief G4Para represents a parallelepiped, essentially a box with half
 * lengths dx,dy,dz 'skewed' so that there are angles theta & phi of the
 * polar line joining the faces at +-dz in z, and alpha formed by the y axis
 * and the plane joining the centre of the faces parallel to the z-x plane
 * at -dy and +dy.
 */

class G4Para : public G4CSGSolid
{
  public:

    /**
     * Constructs a parallelepiped, given a name and its parameters.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx Half-length in x.
     *  @param[in] pDy Half-length in y.
     *  @param[in] pDz Half-length in z.
     *  @param[in] pAlpha Angle formed by the Y axis and by the plane joining
     *             the centre of the faces parallel to the Z-X plane at -dy
     *             and +dy.
     *  @param[in] pTheta Polar angle of the line joining the centres of the
     *             faces at -dz and +dz in Z.
     *  @param[in] pPhi Azimuthal angle of the line joining the centres of
     *             the faces at -dz and +dz in Z.
     */
    G4Para(const G4String& pName,
                 G4double pDx, G4double pDy, G4double pDz,
                 G4double pAlpha, G4double pTheta, G4double pPhi);

    /**
     * Constructs a parallelepiped, given a name and its 8 vertices.
     *  @param[in] pName The name of the solid.
     *  @param[in] pt Points of the 8 vertices.
     */
    G4Para(const G4String& pName,
           const G4ThreeVector pt[8]);

    /**
     * Default destructor.
     */
    ~G4Para() override = default;

    /**
     * Accessors. Obtain (re)computed values of the original parameters.
     */
    inline G4double GetZHalfLength()  const;
    inline G4ThreeVector GetSymAxis() const;
    inline G4double GetYHalfLength()  const;
    inline G4double GetXHalfLength()  const;
    inline G4double GetTanAlpha()     const;
    inline G4double GetAlpha()  const;
    inline G4double GetTheta()  const;    
    inline G4double GetPhi()    const;
   
    /**
     * Modifiers.
     */
    inline void SetXHalfLength(G4double val);
    inline void SetYHalfLength(G4double val);
    inline void SetZHalfLength(G4double val);
    inline void SetAlpha(G4double alpha);
    inline void SetTanAlpha(G4double val);
    inline void SetThetaAndPhi(G4double pTheta, G4double pPhi);
   
    /**
     * Sets all parameters, as for constructor.
     */
    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                          G4double pAlpha, G4double pTheta, G4double pPhi);

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

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
                                 G4double& pMin, G4double& pMax) const override;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    /**
     * Returns the type ID, "G4Para" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    G4bool IsFaceted() const override;

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
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron () const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4Para(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Para(const G4Para& rhs);
    G4Para& operator=(const G4Para& rhs);

  private:

    /**
     * Checks the dimension parameters given in input.
     */
    void CheckParameters();

    /**
     * Sets the side planes.
     */
    void MakePlanes();

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;
   
  private:

    G4double halfCarTolerance;
    G4double fDx,fDy,fDz;
    G4double fTalpha,fTthetaCphi,fTthetaSphi;
    struct { G4double a,b,c,d; } fPlanes[4];
};

#include "G4Para.icc"

#endif

#endif
