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
// $Id: G4Trd.hh,v 1.8 2002-10-28 11:43:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4Trd
//
// Class description:
//
//   A G4Trd is a trapezoid with the x and y dimensions varying along z
//   functions:
//
//   Member Data:
//
//     fDx1    Half-length along x at the surface positioned at -dz
//     fDx2    Half-length along x at the surface positioned at +dz
//     fDy1    Half-length along y at the surface positioned at -dz
//     fDy2    Half-length along y at the surface positioned at +dz
//     fDz     Half-length along z axis

// History:
// 12.01.95 P.Kent: Old prototype code converted to thick geometry
// 17.02.95 P.Kent: Exiting normal return
// 19.08.96 P.Kent, V.Grichine: Fs in accordance with G4Box
// 21.04.97 J.Apostolakis: Added Set Methods
// 19.11.99 V.Grichine: kUndefined was added to Eside enum 
// --------------------------------------------------------------------

#ifndef G4TRD_HH
#define G4TRD_HH

#include "G4CSGSolid.hh"

class G4Trd : public G4CSGSolid 
{
  public:  // with description

    G4Trd( const G4String& pName,
                 G4double pdx1, G4double pdx2,
                 G4double pdy1, G4double pdy2,
                 G4double pdz );
      //
      // Constructs a trapezoid with name, and half lengths

    virtual ~G4Trd();
      //
      // Destructor

    // Accessors

    inline G4double GetXHalfLength1() const;
    inline G4double GetXHalfLength2() const;
    inline G4double GetYHalfLength1() const;
    inline G4double GetYHalfLength2() const;
    inline G4double GetZHalfLength()  const;

    // Modifiers

    inline void SetXHalfLength1(G4double val);
    inline void SetXHalfLength2(G4double val);
    inline void SetYHalfLength1(G4double val);
    inline void SetYHalfLength2(G4double val);
    inline void SetZHalfLength(G4double val);

    // Methods of solid

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax ) const;

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

    void CheckAndSetAllParameters ( G4double pdx1, G4double pdx2,
                                    G4double pdy1, G4double pdy2,
                                    G4double pdz );

    void SetAllParameters ( G4double pdx1, G4double pdx2,
                            G4double pdy1, G4double pdy2,
                            G4double pdz );

    G4GeometryType GetEntityType() const;

    G4std::ostream& StreamInfo( G4std::ostream& os ) const;

    // Visualisation functions

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4NURBS*      CreateNURBS        () const;

  protected:  // without description

    G4ThreeVectorList*
    CreateRotatedVertices( const G4AffineTransform& pTransform ) const;
      //
      // Creates the List of transformed vertices in the format required
      // for G4CSGSolid:: ClipCrossSection and ClipBetweenSections

    G4double fDx1,fDx2,fDy1,fDy2,fDz;

    // Codes for faces (kPX=plus x face,kMY= minus y face etc)

    enum ESide {kUndefined, kPX,kMX,kPY,kMY,kPZ,kMZ};

};

#include "G4Trd.icc"

#endif
