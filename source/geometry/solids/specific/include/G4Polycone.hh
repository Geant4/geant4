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
// G4Polycone
//
// Class description:
//
//   Class implementing a CSG-like type "PCON" Geant 3.21 volume,
//   inherited from  class G4VCSGfaceted:
//
//   G4Polycone( const G4String& name,
//               G4double phiStart,     // initial phi starting angle
//               G4double phiTotal,     // total phi angle
//               G4int numZPlanes,      // number of z planes
//               const G4double zPlane[],  // position of z planes
//               const G4double rInner[],  // tangent distance to inner surface
//               const G4double rOuter[])  // tangent distance to outer surface
//
//   Alternative constructor, but limited to increasing-only Z sections:
//
//   G4Polycone( const G4String& name,
//               G4double phiStart,   // initial phi starting angle
//               G4double phiTotal,   // total phi angle
//               G4int    numRZ,      // number corners in r,z space
//               const G4double r[],  // r coordinate of these corners
//               const G4double z[])  // z coordinate of these corners

// Author: David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------
#ifndef G4POLYCONE_HH
#define G4POLYCONE_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPOLYCONE 1
#endif

#if defined(G4GEOM_USE_UPOLYCONE)
  #define G4UPolycone G4Polycone
  #include "G4UPolycone.hh"
#else

#include "G4VCSGfaceted.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyconeHistorical.hh"
#include "G4Polyhedron.hh"

class G4EnclosingCylinder;
class G4ReduciblePolygon;
class G4VCSGface;

class G4Polycone : public G4VCSGfaceted
{
  public:

    G4Polycone( const G4String& name,
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int numZPlanes,     // number of z planes
                const G4double zPlane[],    // position of z planes
                const G4double rInner[],    // tangent distance to inner surface
                const G4double rOuter[]  ); // tangent distance to outer surface

    G4Polycone( const G4String& name,
                      G4double phiStart,    // initial phi starting angle
                      G4double phiTotal,    // total phi angle
                      G4int    numRZ,       // number corners in r,z space
                const G4double r[],         // r coordinate of these corners
                const G4double z[]       ); // z coordinate of these corners

    ~G4Polycone() override;

    EInside Inside( const G4ThreeVector& p ) const override;
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn( const G4ThreeVector& p ) const override;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const override;

    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    G4ThreeVector GetPointOnSurface() const override;

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    G4GeometryType GetEntityType() const override;

    G4VSolid* Clone() const override;

    std::ostream& StreamInfo(std::ostream& os) const override;

    G4Polyhedron* CreatePolyhedron() const override;

    G4bool Reset();

    // Accessors

    inline G4double GetStartPhi()    const;
    inline G4double GetEndPhi()      const;
    inline G4double GetSinStartPhi() const;
    inline G4double GetCosStartPhi() const;
    inline G4double GetSinEndPhi()   const;
    inline G4double GetCosEndPhi()   const;
    inline G4bool IsOpen()           const;
    inline G4int  GetNumRZCorner()   const;
    inline G4PolyconeSideRZ GetCorner(G4int index) const;
    inline G4PolyconeHistorical* GetOriginalParameters() const;
    inline void SetOriginalParameters(G4PolyconeHistorical* pars);

    G4Polycone(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Polycone( const G4Polycone& source );
    G4Polycone &operator=( const G4Polycone& source );
      // Copy constructor and assignment operator.

  protected:

    // Generic initializer, called by all constructors

    G4bool SetOriginalParameters(G4ReduciblePolygon* rz);

    void Create( G4double phiStart,        // initial phi starting angle
                 G4double phiTotal,        // total phi angle
                 G4ReduciblePolygon* rz ); // r/z coordinate of these corners

    void CopyStuff( const G4Polycone& source );

    // Methods for random point generation

    void SetSurfaceElements() const;

  protected:

    // Here are our parameters

    G4double startPhi;        // Starting phi value (0 < phiStart < 2pi)
    G4double endPhi;          // End phi value (0 < endPhi-phiStart < 2pi)
    G4bool phiIsOpen = false; // True if there is a phi segment
    G4int numCorner;          // Number RZ points
    G4PolyconeSideRZ* corners = nullptr; // Corner r,z points
    G4PolyconeHistorical* original_parameters = nullptr; // Original input pars

    G4EnclosingCylinder* enclosingCylinder = nullptr; // Our quick test

    struct surface_element { G4double area = 0.; G4int i0 = 0, i1 = 0, i2 = 0; };
    mutable std::vector<surface_element>* fElements = nullptr;
};

#include "G4Polycone.icc"

#endif

#endif
