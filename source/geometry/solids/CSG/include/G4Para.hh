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
// $Id: G4Para.hh 104452 2017-05-31 15:41:24Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
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

// History:
// 21.3.94 P.Kent Old C++ code converted to tolerant geometry
// 31.10.96 V.Grichine Modifications according G4Box/Tubs before to commit
// 18.11.99 V.Grichine , kUndefined was added to ESide
// --------------------------------------------------------------------

#ifndef G4Para_HH
#define G4Para_HH

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

    virtual ~G4Para();

    // Accessors

    inline G4double GetZHalfLength()  const;
    inline G4ThreeVector GetSymAxis() const;
    inline G4double GetYHalfLength()  const;
    inline G4double GetXHalfLength()  const;
    inline G4double GetTanAlpha()     const;

    // Modifiers

    inline void SetXHalfLength(G4double val);
    inline void SetYHalfLength(G4double val);
    inline void SetZHalfLength(G4double val);
    inline void SetAlpha(G4double alpha);
    inline void SetTanAlpha(G4double val);
    inline void SetThetaAndPhi(double pTheta, double pPhi);

    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                          G4double pAlpha, G4double pTheta, G4double pPhi);

    // Methods of solid

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4Polyhedron* CreatePolyhedron   () const;

  public:  // without description

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

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p) const;
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
