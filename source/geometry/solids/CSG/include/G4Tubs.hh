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
// $Id: G4Tubs.hh,v 1.7 2002-10-28 11:43:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//   the z axis, about which it is centred, and a given
//   minimum and maximum radius. A minimum radius of 0
//   signifies a filled tube /cylinder. The tube segment is
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

// History:
// 10.08.95 P.Kent: General cleanup, use G4VSolid extent helper functions
//                  to CalculateExtent()
// 23.01.94 P.Kent: Converted to `tolerant' geometry
// 19.07.96 J.Allison: G4GraphicsScene - see G4Box
// 22.07.96 J.Allison: Changed SendPolyhedronTo to CreatePolyhedron
// --------------------------------------------------------------------

#ifndef G4TUBS_HH
#define G4TUBS_HH

#include "G4CSGSolid.hh"

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

    // Modifiers

    inline void SetInnerRadius   (G4double newRMin);
    inline void SetOuterRadius   (G4double newRMax);
    inline void SetZHalfLength   (G4double newDz);
    inline void SetStartPhiAngle (G4double newSPhi);
    inline void SetDeltaPhiAngle (G4double newDPhi);

    // Methods for solid

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

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

    G4std::ostream& StreamInfo( G4std::ostream& os ) const;

    // Visualisation functions

    void                DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron*       CreatePolyhedron   () const;
    G4NURBS*            CreateNURBS        () const;

  public:  // without description

    //  Older names for access functions

    inline G4double GetRMin() const;
    inline G4double GetRMax() const;
    inline G4double GetDz  () const;
    inline G4double GetSPhi() const;
    inline G4double GetDPhi() const;

  protected:

    G4ThreeVectorList*
    CreateRotatedVertices( const G4AffineTransform& pTransform ) const;
      //
      // Creates the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections

    G4double fRMin,fRMax,fDz,fSPhi,fDPhi;

    // Used by distanceToOut

    enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

    // used by normal

    enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

};

#include "G4Tubs.icc"

#endif
