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
// G4Trap
//
// Class description:
//
// A G4Trap is a general trapezoid: The faces perpendicular to the
// z planes are trapezia, and their centres are not necessarily on
// a line parallel to the z axis.
//
// Note that of the 11 parameters described below, only 9 are really
// independent - a check for planarity is made in the calculation of the
// equation for each plane. If the planes are not parallel, a call to
// G4Exception is made.
//
//      pDz     Half-length along the z-axis
//      pTheta  Polar angle of the line joining the centres of the faces
//              at -/+pDz
//      pPhi    Azimuthal angle of the line joining the centre of the face at
//              -pDz to the centre of the face at +pDz
//      pDy1    Half-length along y of the face at -pDz
//      pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
//      pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
//      pAlp1   Angle with respect to the y axis from the centre of the side
//              at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz
//
//      pDy2    Half-length along y of the face at +pDz
//      pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
//      pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
//      pAlp2   Angle with respect to the y axis from the centre of the side
//              at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz
//
//
// Member Data:
//
//      fDz     Half-length along the z axis
//      fTthetaCphi = std::tan(pTheta)*std::cos(pPhi)
//      fTthetaSphi = std::tan(pTheta)*std::sin(pPhi)
//      These combinations are suitable for creation of the trapezoid corners
//
//      fDy1    Half-length along y of the face at -fDz
//      fDx1    Half-length along x of the side at y=-fDy1 of the face at -fDz
//      fDx2    Half-length along x of the side at y=+fDy1 of the face at -fDz
//      fTalpha1   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy1 to the centre at y=+fDy1 of the face
//                 at -fDz
//
//      fDy2    Half-length along y of the face at +fDz
//      fDx3    Half-length along x of the side at y=-fDy2 of the face at +fDz
//      fDx4    Half-length along x of the side at y=+fDy2 of the face at +fDz
//      fTalpha2   Tan of Angle with respect to the y axis from the centre of
//                 the side at y=-fDy2 to the centre at y=+fDy2 of the face
//                 at +fDz
//
//      TrapSidePlane fPlanes[4]   Plane equations of the faces not at +/-fDz
//                                 NOTE: order is important !!!

// Author: Paul Kent, 23.03.1994 - Code converted to tolerant geometry
// --------------------------------------------------------------------
#ifndef G4TRAP_HH
#define G4TRAP_HH

#include "G4Types.hh"

struct TrapSidePlane
{
    G4double a,b,c,d;    // Normal unit vector (a,b,c)  and offset (d)
        // => Ax+By+Cz+D=0
};

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTRAP 1
#endif

#if defined(G4GEOM_USE_UTRAP)
  #define G4UTrap G4Trap
  #include "G4UTrap.hh"
#else

#include "G4CSGSolid.hh"

/**
 * @brief G4Trap is a general trapezoid: the faces perpendicular to the Z
 * planes are trapezia, and their centres are not necessarily on a line parallel
 * to the Z axis. A check for planarity is made in the calculation of the
 * equation for each plane. If the planes are not parallel, a call to
 * G4Exception is made.
 */

class G4Trap : public G4CSGSolid
{
  public:

    /**
     * The most general constructor for G4Trap which prepares plane
     * equations and corner coordinates from parameters.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDz Half-length along the Z-axis.
     *  @param[in] pTheta Polar angle of the line joining the centres
     *             of the faces at -/+pDz.
     *  @param[in] pPhi Azimuthal angle of the line joining the centre
     *             of the face at -pDz to the centre of the face at +pDz.
     *  @param[in] pDy1 Half-length along Y of the face at -pDz.
     *  @param[in] pDx1 Half-length along X of the side at y=-pDy1
     *             of the face at -pDz.
     *  @param[in] pDx2 Half-length along X of the side at y=+pDy1
     *             of the face at -pDz.
     *  @param[in] pAlp1 Angle with respect to the Y axis from the centre of the
     *             side at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz.
     *  @param[in] pDy2 Half-length along Y of the face at +pDz.
     *  @param[in] pDx3 Half-length along X of the side at y=-pDy2
     *             of the face at +pDz.
     *  @param[in] pDx4 Half-length along X of the side at y=+pDy2
     *             of the face at +pDz.
     *  @param[in] pAlp2 Angle with respect to the Y axis from the centre of the
     *             side at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz.
     */
    G4Trap( const G4String& pName,
                  G4double pDz,
                  G4double pTheta, G4double pPhi,
                  G4double pDy1, G4double pDx1, G4double pDx2,
                  G4double pAlp1,
                  G4double pDy2, G4double pDx3, G4double pDx4,
                  G4double pAlp2 );

    /**
     * Prepares plane equations and parameters from corner coordinates.
     *  @param[in] pName The name of the solid.
     *  @param[in] pt Points of the 8 vertices.
     */
    G4Trap( const G4String& pName,
            const G4ThreeVector pt[8] ) ;

    /**
     * Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX).
     *  @param[in] pName The name of the solid.
     *  @param[in] pZ Length along Z.
     *  @param[in] pY Length along Y.
     *  @param[in] pX Length along X at the wider side.
     *  @param[in] pLTX Length along X at the narrower side (plTX<=pX).
     */
    G4Trap( const G4String& pName,
                  G4double pZ,
                  G4double pY,
                  G4double pX, G4double pLTX );

    /**
     * Constructor for G4Trd.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx1 Half-length along X at the surface positioned at -dz.
     *  @param[in] pDx2 Half-length along X at the surface positioned at +dz.
     *  @param[in] pDy1 Half-length along Y at the surface positioned at -dz.
     *  @param[in] pDy2 Half-length along Y at the surface positioned at +dz.
     *  @param[in] pDz Half-length along Z axis.
     */
    G4Trap( const G4String& pName,
                  G4double pDx1,  G4double pDx2,
                  G4double pDy1,  G4double pDy2,
                  G4double pDz );

    /**
     * Constructor for G4Para.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx Half-length in X.
     *  @param[in] pDy Half-length in Y.
     *  @param[in] pDz Half-length in Z.
     *  @param[in] pAlpha Angle formed by the Y axis and the plane joining the
     *             centre of the faces parallel to the Z-X plane at -dy and +dy.
     *  @param[in] pTheta Polar angle of the line joining the centres of the
     *             faces at -dz and +dz in Z.
     *  @param[in] pPhi Azimuthal angle of the line joining the centres of
     *             the faces at -dz and +dz in Z.
     */
     G4Trap(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi );

    /**
     * Constructor for "nominal" G4Trap whose parameters are to be set
     * by a G4VPVParamaterisation later on.
     *  @param[in] pName The name of the solid.
     */
     G4Trap( const G4String& pName );

    /**
     * Default destructor.
     */
     ~G4Trap() override = default;

    /**
     * Accessors. Returning the coordinates of a unit vector along a straight
     * line joining centers of -/+fDz planes.
     */
    inline G4double GetZHalfLength()  const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetTanAlpha1()    const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetXHalfLength3() const;
    inline G4double GetXHalfLength4() const;
    inline G4double GetTanAlpha2()    const;

    /**
     * More accessors.
     */
    inline TrapSidePlane GetSidePlane( G4int n ) const;
    inline G4ThreeVector GetSymAxis() const;

    /**
     * Accessors obtaining (re)computed values of the original parameters.
     */
    inline G4double GetPhi() const;
    inline G4double GetTheta() const;
    inline G4double GetAlpha1() const;
    inline G4double GetAlpha2() const;   
   
    /**
     * Sets all parameters, as for constructor. Checks and sets half-widths
     * as well as angles. Makes a final check of co-planarity.
     */
    void SetAllParameters ( G4double pDz,
                            G4double pTheta,
                            G4double pPhi,
                            G4double pDy1,
                            G4double pDx1,
                            G4double pDx2,
                            G4double pAlp1,
                            G4double pDy2,
                            G4double pDx3,
                            G4double pDx4,
                            G4double pAlp2 );

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
    void ComputeDimensions( G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

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
    EInside Inside( const G4ThreeVector& p ) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn( const G4ThreeVector& p ) const override;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    /**
     * Returns the type ID, "G4Trap" of the solid.
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
    std::ostream& StreamInfo( std::ostream& os ) const override;

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
    G4Trap(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4Trap(const G4Trap& rhs);
    G4Trap& operator=(const G4Trap& rhs);

  protected:

    /**
     * Internal methods for checking and building planes.
     * Computing the vertices and setting side planes, checking for planarity.
     */
    void MakePlanes();
    void MakePlanes( const G4ThreeVector pt[8] );

    /**
     * Calculates the coefficents of the plane p1->p2->p3->p4->p1
     * where the ThreeVectors 1-4 are in anti-clockwise order when viewed
     * from infront of the plane (i.e. from normal direction).
     *  @return true if the points are co-planar, false otherwise.
     */
    G4bool MakePlane( const G4ThreeVector& p1,
                      const G4ThreeVector& p2,
                      const G4ThreeVector& p3,
                      const G4ThreeVector& p4,
                            TrapSidePlane& plane ) ;
    /**
     * Recomputes parameters using planes.
     */
    void SetCachedValues();

  private:

    /**
     * Checks the input parameters.
     */
    void CheckParameters();

    /**
     * Computes the coordinates of the trap vertices from planes.
     */
    void GetVertices(G4ThreeVector pt[8]) const;

    /**
     * Algorithm for SurfaceNormal() following the original specification
     * for points not on the surface.
     */
    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;

  private:

    G4double halfCarTolerance;
    G4double fDz,fTthetaCphi,fTthetaSphi;
    G4double fDy1,fDx1,fDx2,fTalpha1;
    G4double fDy2,fDx3,fDx4,fTalpha2;
    TrapSidePlane fPlanes[4];
    G4double fAreas[6];
    G4int fTrapType;
};

#include "G4Trap.icc"

#endif

#endif
