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
//   A G4Parallepiped, essentially a box with half lengths dx,dy,dz
//   `skewed' so that there are angles theta & phi of the polar line
//   joining the faces at +-dz in z, and alpha formed by the y axis
//   and the plane joinng the centre of the faces G4Parallel to the
//   z-x plane at -dy and +dy.
//
//   A G4Para is defined by:
//   dx,dy,dz - Half-length in x,y,z
//   alpha    - Angle formed by the y axis and by the plane joining
//              the centre of the faces G4Parallel to the z-x plane
//              at -dy and +dy
//   theta    - Polar angle of the line joining the centres of the
//              faces at -dz and +dz in z
//   phi      - Azimuthal angle of the line joining the centres of the
//              faces at -dz and +dz in z
//   Member data:
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

// 21.3.94 P.Kent Old C++ code converted to tolerant geometry
// 31.10.96 V.Grichine Modifications according G4Box/Tubs before to commit
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

class G4Para : public G4CSGSolid
{
  public:  // with description

    G4Para(const G4String& pName,
                 G4double pDx, G4double pDy, G4double pDz,
                 G4double pAlpha, G4double pTheta, G4double pPhi);

    G4Para(const G4String& pName,
           const G4ThreeVector pt[8]);

    ~G4Para() override;

    // Accessors

    inline G4double GetZHalfLength()  const;
    inline G4ThreeVector GetSymAxis() const;
    inline G4double GetYHalfLength()  const;
    inline G4double GetXHalfLength()  const;
    inline G4double GetTanAlpha()     const;

    inline G4double GetAlpha()  const;
    inline G4double GetTheta()  const;    
    inline G4double GetPhi()    const;
    // Obtain (re)computed values of original parameters
   
    // Modifiers

    inline void SetXHalfLength(G4double val);
    inline void SetYHalfLength(G4double val);
    inline void SetZHalfLength(G4double val);
    inline void SetAlpha(G4double alpha);
    inline void SetTanAlpha(G4double val);
    inline void SetThetaAndPhi(G4double pTheta, G4double pPhi);
   
    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                          G4double pAlpha, G4double pTheta, G4double pPhi);

    // Methods of solid

    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const override;

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

    G4GeometryType GetEntityType() const override;

    G4ThreeVector GetPointOnSurface() const override;

    G4VSolid* Clone() const override;

    std::ostream& StreamInfo(std::ostream& os) const override;

    // Visualisation functions

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron   () const override;

    G4Para(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects

    G4Para(const G4Para& rhs);
    G4Para& operator=(const G4Para& rhs);
      // Copy constructor and assignment operator

  private:

    void CheckParameters();
      // Check parameters

    void MakePlanes();
      // Set side planes

    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

   
  private:

    G4double halfCarTolerance;
    G4double fDx,fDy,fDz;
    G4double fTalpha,fTthetaCphi,fTthetaSphi;
    struct { G4double a,b,c,d; } fPlanes[4];
};

#include "G4Para.icc"

#endif

#endif
