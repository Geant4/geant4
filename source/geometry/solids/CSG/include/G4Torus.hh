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
// $Id: G4Torus.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Torus
//
// Class description:
//
//   A torus or torus segment with curved sides parallel to the z-axis.
//   The torus has a specified swept radius about which it is centered,
//   and a given minimum and maximum radius. A minimum radius of 0
//   signifies a filled torus.
//   The torus segment is specified by starting and delta angles for phi,
//   with 0 being the +x axis, PI/2 the +y axis. A delta angle of 2PI
//   signifies a complete, unsegmented torus/cylindr.
//
//   Member functions:
//
//   As inherited from G4CSGSolid+
//
//     G4Torus(const G4String      &pName
//             G4double      pRmin
//             G4double      pRmax
//             G4double      pRtor
//             G4double      pSPhi
//             G4double      pDPhi )
//
//     - Construct a torus with the given name and dimensions.
//       The angles are provided is radians. pRtor >= pRmax
//
//   Member Data:
//
//  fRmin  Inside radius
//  fRmax  Outside radius
//  fRtor  swept radius of torus
//
//  fSPhi  The starting phi angle in radians,
//         adjusted such that fSPhi+fDPhi<=2PI, fSPhi>-2PI
//
//  fDPhi  Delta angle of the segment in radians
//
//   You could find very often in G4Torus functions values like 'pt' or
//   'it'. These are the distances from p or i G4ThreeVector points in the
//   plane (Z axis points p or i) to fRtor point in XY plane. This value is
//   similar to rho for G4Tubs and is used for definiton of the point
//   relative to fRmin and fRmax, i.e. for solution of inside/outside
//   problems
 
// History:
// 30.10.96 V.Grichine: first version of G4Torus
// 21.04.98 J.Apostolakis: added SetAllParameters() function
// 26.05.00 V.Grichine: added new SolveBiQuadratic/Cubic() developed
//                      by O.Cremonesi 
// 31.08.00 E.Medernach: added SolveNumeric functions, migrated to
//                       numeric solutions
// --------------------------------------------------------------------

#ifndef G4TORUS_HH
#define G4TORUS_HH

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTORUS 1
#endif

#if (defined(G4GEOM_USE_UTORUS) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UTorus G4Torus
  #include "G4UTorus.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"

class G4Torus : public G4CSGSolid
{

  public:  // with description

    G4Torus(const G4String &pName,
                  G4double pRmin,
                  G4double pRmax,
                  G4double pRtor,
                  G4double pSPhi,
                  G4double pDPhi);

   ~G4Torus();
    
    // Accessors

    inline G4double GetRmin() const;
    inline G4double GetRmax() const;
    inline G4double GetRtor() const;
    inline G4double GetSPhi() const;
    inline G4double GetDPhi() const;
    inline G4double GetSinStartPhi () const;
    inline G4double GetCosStartPhi () const;
    inline G4double GetSinEndPhi   () const;
    inline G4double GetCosEndPhi   () const;

    // Methods of solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    EInside Inside(const G4ThreeVector& p) const;
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;
    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;
    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions

    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4Polyhedron*       CreatePolyhedron   () const;

  public:  // without description

    void SetAllParameters(G4double pRmin, G4double pRmax, G4double pRtor,
                          G4double pSPhi, G4double pDPhi);
 
    G4Torus(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Torus(const G4Torus& rhs);
    G4Torus& operator=(const G4Torus& rhs); 
      // Copy constructor and assignment operator.

  private:

    void TorusRootsJT(const G4ThreeVector& p,
                      const G4ThreeVector& v,
                            G4double r,
                            std::vector<G4double>& roots) const ;

    G4double SolveNumericJT(const G4ThreeVector& p,
                            const G4ThreeVector& v,
                                  G4double r,
                                  G4bool IsDistanceToIn) const;

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p) const;
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

  private:

    G4double fRmin,fRmax,fRtor,fSPhi,fDPhi;

    // Used by distanceToOut
    enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi};

    // used by normal
    enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi};
    
    G4double fRminTolerance, fRmaxTolerance, kRadTolerance, kAngTolerance;
      // Radial and angular tolerances

    G4double halfCarTolerance, halfAngTolerance;
      // Cached half tolerance values

};

#include "G4Torus.icc"

#endif  // defined(G4GEOM_USE_UTORUS) && defined(G4GEOM_USE_SYS_USOLIDS)


#endif // G4TORUS_HH
