//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4EllipticalTube.hh,v 1.13 2004/10/10 10:39:11 johna Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

class G4EllipticalTube : public G4VSolid
{
  public:  // with description

    G4EllipticalTube( const G4String &name, 
                            G4double theDx,
                            G4double theDy,
                            G4double theDz );

    virtual ~G4EllipticalTube();

    // Standard solid methods

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

    std::ostream& StreamInfo(std::ostream& os) const;

    G4double GetCubicVolume();

    // Visualisation methods

    G4Polyhedron* CreatePolyhedron() const;
    void DescribeYourselfTo( G4VGraphicsScene& scene ) const;
    G4VisExtent GetExtent() const;
    virtual G4Polyhedron* GetPolyhedron () const;

    // Accessors

    inline G4double GetDx() const;
    inline G4double GetDy() const;
    inline G4double GetDz() const;
  
    inline void SetDx( const G4double newDx );
    inline void SetDy( const G4double newDy );
    inline void SetDz( const G4double newDz );
 
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

    G4double fCubicVolume;
    mutable G4Polyhedron* fpPolyhedron;

};

#include "G4EllipticalTube.icc"

#endif
