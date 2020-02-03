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
//   A G4Trap is a general trapezoid: The faces perpendicular to the
//   z planes are trapezia, and their centres are not necessarily on
//   a line parallel to the z axis.
//
//   Note that of the 11 parameters described below, only 9 are really
//   independent - a check for planarity is made in the calculation of the
//   equation for each plane. If the planes are not parallel, a call to
//   G4Exception is made.
//
//      pDz     Half-length along the z-axis
//      pTheta  Polar angle of the line joining the centres of the faces
//              at -/+pDz
//      pPhi    Azimuthal angle of the line joing the centre of the face at
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
//   Member Data:
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

// 23.3.94 P.Kent: Old C++ code converted to tolerant geometry
// 9.9.96  V.Grichine: Final modifications before to commit
// 8.12.97 J.Allison: Added "nominal" contructor and method SetAllParameters
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

class G4Trap : public G4CSGSolid
{

  public:  // with description

    G4Trap( const G4String& pName,
                  G4double pDz,
                  G4double pTheta, G4double pPhi,
                  G4double pDy1, G4double pDx1, G4double pDx2,
                  G4double pAlp1,
                  G4double pDy2, G4double pDx3, G4double pDx4,
                  G4double pAlp2 );
      //
      // The most general constructor for G4Trap which prepares plane
      // equations and corner coordinates from parameters

    G4Trap( const G4String& pName,
            const G4ThreeVector pt[8] ) ;
      //
      // Prepares plane equations and parameters from corner coordinates

    G4Trap( const G4String& pName,
                  G4double pZ,
                  G4double pY,
                  G4double pX, G4double pLTX );
      //
      // Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX)

    G4Trap( const G4String& pName,
                  G4double pDx1,  G4double pDx2,
                  G4double pDy1,  G4double pDy2,
                  G4double pDz );
      //
      // Constructor for G4Trd

     G4Trap(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi );
      //
      // Constructor for G4Para

     G4Trap( const G4String& pName );
       //
       // Constructor for "nominal" G4Trap whose parameters are to be set
       // by a G4VPVParamaterisation later

     virtual ~G4Trap() ;
       //
       // Destructor

  // Accessors

    inline G4double GetZHalfLength()  const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetTanAlpha1()    const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetXHalfLength3() const;
    inline G4double GetXHalfLength4() const;
    inline G4double GetTanAlpha2()    const;
      //
      // Returns coordinates of unit vector along straight
      // line joining centers of -/+fDz planes

    inline TrapSidePlane GetSidePlane( G4int n ) const;
    inline G4ThreeVector GetSymAxis() const;

  // Modifiers

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

  // Methods for solid

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax ) const;

    EInside Inside( const G4ThreeVector& p ) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;

    G4double DistanceToIn( const G4ThreeVector& p ) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const;

    G4double DistanceToOut( const G4ThreeVector& p ) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo( std::ostream& os ) const;

  // Visualisation functions

    void          DescribeYourselfTo ( G4VGraphicsScene& scene  ) const;
    G4Polyhedron* CreatePolyhedron   () const;

  public:  // without description

    G4Trap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Trap(const G4Trap& rhs);
    G4Trap& operator=(const G4Trap& rhs);
      // Copy constructor and assignment operator.

  protected:  // with description

    void MakePlanes();
    void MakePlanes(const G4ThreeVector pt[8]);
    G4bool MakePlane( const G4ThreeVector& p1,
                      const G4ThreeVector& p2,
                      const G4ThreeVector& p3,
                      const G4ThreeVector& p4,
                            TrapSidePlane& plane ) ;

  private:

    void CheckParameters();
      // Check parameters

    void GetVertices(G4ThreeVector pt[8]) const;
      // Compute coordinates of the trap vertices from planes

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

  private:

    G4double halfCarTolerance;
    G4double fDz,fTthetaCphi,fTthetaSphi;
    G4double fDy1,fDx1,fDx2,fTalpha1;
    G4double fDy2,fDx3,fDx4,fTalpha2;
    TrapSidePlane fPlanes[4];
    G4int fTrapType;
};

#include "G4Trap.icc"

#endif

#endif
