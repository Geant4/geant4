// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Tubs.hh,v 1.3 2000-04-07 12:55:03 gcosmo Exp $
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
//   Member functions:
//
//   As inherited from G4VSolid+
//
//     G4Tubs(const G4String      &pName
//            G4double      pRMin
//            G4double      pRMax
//            G4double      pDz
//            G4double      pSPhi
//            G4double      pDPhi )
//
//       - Construct a tubs with the given name and dimensions.
//         The angles are provided is radians.
//
//
//   Protected:
//
//     G4ThreeVectorList*
//     CreateRotatedVertices(const G4AffineTransform& pTransform) const
//
//       - Create the List of transformed vertices in the format required
//         for G4VSolid:: ClipCrossSection and ClipBetweenSections.
//   
//   Member Data:
//
//	fRMin	Inner radius
//	fRMax	Outer radius
//	fDz	half length in z
//
//	fSPhi	The starting phi angle in radians,
//              adjusted such the fSPhi+fDPhi<=2PI,
//              fSPhi>-2PI
//
//	fDPhi	Delta angle of the segment in radians

// History:
// 10.8.95 P.Kent General cleanup, use G4VSolid extent helper functions
//                to CaluclateExtent
// 23.1.94 P.Kent Converted to `tolerant' geometry
// 19.07.96 J.Allison G4GraphicsScene - see G4Box.
// 22.07.96 J.Allison Changed SendPolyhedronTo to CreatePolyhedron.
// --------------------------------------------------------------------

#ifndef G4TUBS_HH
#define G4TUBS_HH

#include "G4CSGSolid.hh"

class G4Tubs : public G4CSGSolid {
public:
    G4Tubs(const G4String &pName,
	   G4double pRMin,
	   G4double pRMax,
	   G4double pDz,
	   G4double pSPhi,
	   G4double pDPhi);

    virtual ~G4Tubs();
    
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
			   const G4VoxelLimits& pVoxelLimit,
			   const G4AffineTransform& pTransform,
			   G4double& pmin, G4double& pmax) const;

    G4double    GetInnerRadius  () const { return fRMin; }
    G4double    GetOuterRadius () const { return fRMax; }
    G4double    GetZHalfLength   () const { return fDz  ; }
    G4double    GetStartPhiAngle () const { return fSPhi; }
    G4double    GetDeltaPhiAngle () const { return fDPhi; }

    void        SetInnerRadius  (G4double newRMin) { fRMin= newRMin; }
    void        SetOuterRadius (G4double newRMax) { fRMax= newRMax; } 
    void        SetZHalfLength   (G4double newDz)   { fDz=   newDz  ; }
    void        SetStartPhiAngle (G4double newSPhi) { fSPhi= newSPhi; }
    void        SetDeltaPhiAngle (G4double newDPhi) { fDPhi= newDPhi; }

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			   const G4bool calcNorm=G4bool(false),
			   G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    // Naming method (pseudo-RTTI : run-time type identification)
    virtual G4GeometryType  GetEntityType() const { return G4String("G4Tubs"); }

    // Visualisation functions
    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent         GetExtent          () const;
    G4Polyhedron*       CreatePolyhedron   () const;
    G4NURBS*            CreateNURBS        () const;

    //  Older names for access functions
    G4double    GetRMin() const { return GetInnerRadius(); }   // fRMin 
    G4double    GetRMax() const { return GetOuterRadius(); }  // fRMax 
    G4double    GetDz  () const { return GetZHalfLength()  ; }  // fDz   
    G4double    GetSPhi() const { return GetStartPhiAngle(); }  // fSPhi 
    G4double    GetDPhi() const { return GetDeltaPhiAngle(); }  // fDPhi

protected:
    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;

    G4double fRMin,fRMax,fDz,fSPhi,fDPhi;

  // Used by distanceToOut
  enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};
  // used by normal
  enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

};
   	
#endif
