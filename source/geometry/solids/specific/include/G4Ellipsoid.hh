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
// $Id: G4Ellipsoid.hh,v 1.9 2005/11/09 15:04:28 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Ellipsoid
//
// Class description:
//
//   A G4Ellipsoid is an ellipsoidal solid, optionally cut at a given z.
//
//   Member Data:
//
//      xSemiAxis       semi-axis, x
//      ySemiAxis       semi-axis, y
//      zSemiAxis       semi-axis, z
//      zBottomCut      lower cut plane level, z (solid lies above this plane)
//      zTopCut         upper cut plane level, z (solid lies below this plane)

// History:
// -------
// 10.11.1999  G.Horton-Smith (Caltech, USA) - First implementation
// 10.02.2005  G.Guerrieri (INFN Genova, Italy) - Revision
// --------------------------------------------------------------------
#ifndef G4Ellipsoid_HH
#define G4Ellipsoid_HH

#include "G4VSolid.hh"

class G4Ellipsoid : public G4VSolid
{
  public:  // with description

    G4Ellipsoid(const G4String& pName,
                      G4double  pxSemiAxis,
                      G4double  pySemiAxis,
                      G4double  pzSemiAxis,
                      G4double  pzBottomCut=0,
                      G4double  pzTopCut=0);

    virtual ~G4Ellipsoid();

    // Access functions
   
    inline G4double GetSemiAxisMax (G4int i) const;
    inline G4double GetZBottomCut() const;
    inline G4double GetZTopCut() const;
    inline void SetSemiAxis (G4double x, G4double y, G4double z);
    inline void SetZCuts (G4double newzBottomCut, G4double newzTopCut);
    inline G4double GetCubicVolume(); 

    // Solid standard methods
   
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;
    EInside Inside(const G4ThreeVector& p) const;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm=G4bool(false),
                                 G4bool *validNorm=0,
                                 G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    G4ThreeVector GetPointOnSurface() const;

    // Visualisation functions
  
    G4Polyhedron* GetPolyhedron () const;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent() const;
    G4Polyhedron* CreatePolyhedron() const;
    G4NURBS*      CreateNURBS() const;
       
  public:  // without description

    G4Ellipsoid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  protected:  // without description
 
    G4ThreeVectorList* CreateRotatedVertices(const G4AffineTransform& pT,
                                                   G4int& noPV) const;

    mutable G4Polyhedron* fpPolyhedron;

  private:

    G4double fCubicVolume;
    G4double xSemiAxis, ySemiAxis, zSemiAxis,
             semiAxisMax, zBottomCut, zTopCut;
};

#include "G4Ellipsoid.icc"

#endif
