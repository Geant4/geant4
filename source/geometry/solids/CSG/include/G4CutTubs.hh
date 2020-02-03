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
// G4CutTubs
//
// Class description:
//
// G4CutTubs is a tube with possible cuts in +-Z.
// Implementation adapted from G4Tubs (subclass of G4Tubs) and
// from TGEo Ctube implementation (by A.Gheata, CERN)
//
// G4CutTubs(pName,pRMin,pRMax,pDZ,pSPhi,pEPhi,pLowNorm,pHighNorm)
//           pName,pRMin,pRMax,pDZ,pSPhi,pEPhi are the same as for G4Tubs,
//           pLowNorm=Outside Normal at -Z
//           pHighNorm=Outsie Normal at +Z.

// Author: Tatiana Nikitina, CERN
// --------------------------------------------------------------------

#ifndef G4CUTTUBS_HH
#define G4CUTTUBS_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UCTUBS 1
#endif

#if defined(G4GEOM_USE_UCTUBS)
  #define G4UCutTubs G4CutTubs
  #include "G4UCutTubs.hh"
#else

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

class G4CutTubs : public G4CSGSolid
{
  public:  // with description

    G4CutTubs( const G4String& pName,
                     G4double pRMin,
                     G4double pRMax,
                     G4double pDz,
                     G4double pSPhi,
                     G4double pDPhi,
                     G4ThreeVector pLowNorm,
                     G4ThreeVector pHighNorm );
      //
      // Constructs a tubs with the given name and dimensions

   ~G4CutTubs();
      //
      // Destructor

    // Accessors

    inline G4double GetInnerRadius   () const;
    inline G4double GetOuterRadius   () const;
    inline G4double GetZHalfLength   () const;
    inline G4double GetStartPhiAngle () const;
    inline G4double GetDeltaPhiAngle () const;
    inline G4double GetSinStartPhi   () const;
    inline G4double GetCosStartPhi   () const;
    inline G4double GetSinEndPhi     () const;
    inline G4double GetCosEndPhi     () const;
    inline G4ThreeVector GetLowNorm  () const;
    inline G4ThreeVector GetHighNorm () const;

    // Modifiers

    inline void SetInnerRadius   (G4double newRMin);
    inline void SetOuterRadius   (G4double newRMax);
    inline void SetZHalfLength   (G4double newDz);
    inline void SetStartPhiAngle (G4double newSPhi, G4bool trig=true);
    inline void SetDeltaPhiAngle (G4double newDPhi);

    // Methods for solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pmin, G4double& pmax ) const;

    EInside Inside( const G4ThreeVector& p ) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo( std::ostream& os ) const;

    // Visualisation functions

    void                DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron*       CreatePolyhedron   () const;

  public:  // without description

    G4CutTubs(__void__&);
      //
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4CutTubs(const G4CutTubs& rhs);
    G4CutTubs& operator=(const G4CutTubs& rhs);
      // Copy constructor and assignment operator.

    //  Older names for access functions

    inline G4double GetRMin() const;
    inline G4double GetRMax() const;
    inline G4double GetDz  () const;
    inline G4double GetSPhi() const;
    inline G4double GetDPhi() const;

  protected:

    inline void Initialize();
      //
      // Reset relevant values to zero

    inline void CheckSPhiAngle(G4double sPhi);
    inline void CheckDPhiAngle(G4double dPhi);
    inline void CheckPhiAngles(G4double sPhi, G4double dPhi);
      //
      // Reset relevant flags and angle values

    inline void InitializeTrigonometry();
      //
      // Recompute relevant trigonometric values and cache them

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;
      //
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

    G4bool IsCrossingCutPlanes() const;
      // Check if the cutted planes are crossing.
      // If 'true' , solid is ill defined

    G4double GetCutZ(const G4ThreeVector& p) const;
      // Get Z value of the point on Cutted Plane

    void GetMaxMinZ(G4double& zmin,G4double& zmax)const;
      // Get Max and Min values of Z on Cutted Plane,
      // Used for Calculate BoundingLimits()

  private:

    G4double kRadTolerance, kAngTolerance;
      //
      // Radial and angular tolerances

    G4double fRMin, fRMax, fDz, fSPhi, fDPhi;
      //
      // Radial and angular dimensions

    G4double sinCPhi, cosCPhi, cosHDPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;
      //
      // Cached trigonometric values

    G4bool fPhiFullCutTube = false;
      //
      // Flag for identification of section or full tube

    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;
      //
      // Cached half tolerance values

    G4ThreeVector fLowNorm, fHighNorm;
      //
      // Normals of Cut at -/+ Dz
};

#include "G4CutTubs.icc"

#endif

#endif
