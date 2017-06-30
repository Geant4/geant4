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
// $Id: G4EllipticalTube.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4EllipticalTube
//
// Class description:
//
//   Declaration of a CSG volume representing a tube with elliptical
//   cross section (geant3 solid 'ELTU'):
//   
//   G4EllipticalTube( const G4String& name, 
//                           G4double  Dx,
//                           G4double  Dy,
//                           G4double  Dz )
//
//   The equation of the surface in x/y is 1.0 = (x/dx)**2 + (y/dy)**2

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4EllipticalTube_hh
#define G4EllipticalTube_hh

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

class G4EllipticalTube : public G4VSolid
{
  public:  // with description

    G4EllipticalTube( const G4String &name, 
                            G4double theDx,
                            G4double theDy,
                            G4double theDz );

    virtual ~G4EllipticalTube();

    // Standard solid methods

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pmin, G4double& pmax ) const;
  
    EInside Inside( const G4ThreeVector& p ) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const;
    G4double DistanceToIn( const G4ThreeVector& p ) const;
    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm=false,
                                  G4bool *validNorm=0,
                                  G4ThreeVector *n=0 ) const;
    G4double DistanceToOut( const G4ThreeVector& p ) const;

    G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    G4ThreeVector GetPointOnSurface() const;

    // Visualisation methods

    G4Polyhedron* CreatePolyhedron() const;
    G4Polyhedron* GetPolyhedron () const;
    void DescribeYourselfTo( G4VGraphicsScene& scene ) const;
    G4VisExtent GetExtent() const;

    // Accessors

    inline G4double GetDx() const;
    inline G4double GetDy() const;
    inline G4double GetDz() const;
  
    inline void SetDx( const G4double newDx );
    inline void SetDy( const G4double newDy );
    inline void SetDz( const G4double newDz );
 
  public:  // without description

    G4EllipticalTube(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4EllipticalTube(const G4EllipticalTube& rhs);
    G4EllipticalTube& operator=(const G4EllipticalTube& rhs); 
      // Copy constructor and assignment operator.

  protected:  // without description

    G4double dx, dy, dz;
  
    // Utility

    inline G4double CheckXY( const G4double x,
                             const G4double y,
                             const G4double toler ) const;
    inline G4double CheckXY( const G4double x, const G4double y ) const;

    G4int IntersectXY( const G4ThreeVector &p,
                       const G4ThreeVector &v, G4double s[2] ) const;

  private:

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p ) const;
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface.

    G4double halfTol;

    G4double fCubicVolume;
    G4double fSurfaceArea;
    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;
};

#include "G4EllipticalTube.icc"

#endif
