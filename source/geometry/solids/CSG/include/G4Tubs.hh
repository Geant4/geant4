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
// $Id: G4Tubs.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4Tubs
//
// Class description:
//
//   A tube or tube segment with curved sides parallel to
//   the z-axis. The tube has a specified half-length along
//   the z-axis, about which it is centered, and a given
//   minimum and maximum radius. A minimum radius of 0
//   corresponds to filled tube /cylinder. The tube segment is
//   specified by starting and delta angles for phi, with 0
//   being the +x axis, PI/2 the +y axis.
//   A delta angle of 2PI signifies a complete, unsegmented
//   tube/cylinder.
//
//   Member Data:
//
//   fRMin  Inner radius
//   fRMax  Outer radius
//   fDz  half length in z
//
//   fSPhi  The starting phi angle in radians,
//          adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI
//
//   fDPhi  Delta angle of the segment.
//
//   fPhiFullTube   Boolean variable used for indicate the Phi Section

// History:
// 10.08.95 P.Kent: General cleanup, use G4VSolid extent helper functions
//                  to CalculateExtent()
// 23.01.94 P.Kent: Converted to `tolerant' geometry
// 19.07.96 J.Allison: G4GraphicsScene - see G4Box
// 22.07.96 J.Allison: Changed SendPolyhedronTo to CreatePolyhedron
// --------------------------------------------------------------------
#ifndef G4TUBS_HH
#define G4TUBS_HH

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTUBS 1
#endif

#if defined(G4GEOM_USE_UTUBS)
  #define G4UTubs G4Tubs
  #include "G4UTubs.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

class G4Tubs : public G4CSGSolid
{
  public:  // with description

    G4Tubs( const G4String& pName,
                  G4double pRMin,
                  G4double pRMax,
                  G4double pDz,
                  G4double pSPhi,
                  G4double pDPhi );
      //
      // Constructs a tubs with the given name and dimensions

    virtual ~G4Tubs();
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

    // Modifiers

    inline void SetInnerRadius   (G4double newRMin);
    inline void SetOuterRadius   (G4double newRMax);
    inline void SetZHalfLength   (G4double newDz);
    inline void SetStartPhiAngle (G4double newSPhi, G4bool trig=true);
    inline void SetDeltaPhiAngle (G4double newDPhi);
    
    // Methods for solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

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
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo( std::ostream& os ) const;

    // Visualisation functions

    void                DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron*       CreatePolyhedron   () const;

  public:  // without description

    G4Tubs(__void__&);
      //
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Tubs(const G4Tubs& rhs);
    G4Tubs& operator=(const G4Tubs& rhs); 
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

    virtual G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;
      //
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

  protected:

    // Used by distanceToOut
    //
    enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

    // Used by normal
    //
    enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

    G4double kRadTolerance, kAngTolerance;
      //
      // Radial and angular tolerances

    G4double fRMin, fRMax, fDz, fSPhi, fDPhi;
      //
      // Radial and angular dimensions

    G4double sinCPhi, cosCPhi, cosHDPhiOT, cosHDPhiIT,
             sinSPhi, cosSPhi, sinEPhi, cosEPhi;
      //
      // Cached trigonometric values

    G4bool fPhiFullTube;
      //
      // Flag for identification of section or full tube

    G4double halfCarTolerance, halfRadTolerance, halfAngTolerance;
      //
      // Cached half tolerance values
};

#include "G4Tubs.icc"

#endif

#endif
