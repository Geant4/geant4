// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trd.hh,v 1.4 2000-04-07 12:55:03 gcosmo Exp $
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
//   A Trd is a trapezoid with the x and y dimensions varying along z
//   functions:
//
//   As inherited from G4CSGSolid +
//
//     G4Trd( const G4String& pName,
//            G4double pdx1, G4double pdx2,
//            G4double pdy1, G4double pdy2,
//            G4double pdz )
//
//       - Construct a trapezoid with name, and half lengths
//         dpx1,dpx2,dpy1,dpy2,dpz
//
//     G4double GetXHalfLength1() const
//     G4double GetXHalfLength2() const
//     G4double GetYHalfLength1() const
//     G4double GetYHalfLength2() const
//     G4double GetZHalfLength()  const
//
//       - Return the respective parameter
//
//     void SetXHalfLength1(G4double) 
//     void SetXHalfLength2(G4double)
//     void SetYHalfLength1(G4double)
//     void SetYHalfLength2(G4double)
//     void SetZHalfLength(G4double)
//
//       - Set the respective parameter
//
//   Protected:
// 
//     G4ThreeVectorList*
//     CreateRotatedVertices(const G4AffineTransform& pTransform) const
// 
//       - Create the List of transformed vertices in the format required
//         for G4CSGSolid:: ClipCrossSection and ClipBetweenSections.
//
//
//   Member Data:
//
//     fDx1    Half-length along x at the surface positioned at -dz
//     fDx2    Half-length along x at the surface positioned at +dz
//     fDy1    Half-length along y at the surface positioned at -dz
//     fDy2    Half-length along y at the surface positioned at +dz
//     fDz     Half-length along z axis

// History:
// 19.11.99 V.Grichine, kUndefined was added to Eside enum 
// 21.04.97 J. Apostolakis         Added Set Methods.
// 19.08.96 P. Kent, V. Grichine ->Fs in accordance with G4Box
// 17.02.95 P.Kent Exiting normal return
// 12.01.95 P.Kent Old prototype code converted to thick geometry
// --------------------------------------------------------------------

#ifndef G4TRD_HH
#define G4TRD_HH

#include "G4CSGSolid.hh"

class G4Trd : public G4CSGSolid 
{
  public:

    G4Trd( const G4String& pName,
           G4double pdx1, G4double pdx2,
	   G4double pdy1, G4double pdy2,
           G4double pdz);

    virtual ~G4Trd();

// Access functions
    G4double GetXHalfLength1() const
    {
	return fDx1;
    }

    G4double GetXHalfLength2() const
    {
	return fDx2;
    }
    
    G4double GetYHalfLength1() const
    {
	return fDy1;
    }

    G4double GetYHalfLength2() const
    {
	return fDy2;
    }
    
    G4double GetZHalfLength() const
    {
	return fDz;
    }
    
    void SetXHalfLength1(G4double val) 
    {
	fDx1= val;
    }
    void SetXHalfLength2(G4double val)
    {
	fDx2= val;
    }
    void SetYHalfLength1(G4double val) 
    {
	fDy1= val;
    }
    void SetYHalfLength2(G4double val)
    {
	fDy2= val;
    }
    void SetZHalfLength(G4double val)
    {
	fDz= val;
    }

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
			   const G4VoxelLimits& pVoxelLimit,
			   const G4AffineTransform& pTransform,
			   G4double& pMin, G4double& pMax) const;


    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;

    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
			   const G4bool calcNorm=false,
			   G4bool *validNorm=0,
                           G4ThreeVector *n=0) const;

    G4double DistanceToOut(const G4ThreeVector& p) const;

             // Naming method (pseudo-RTTI : run-time type identification

    virtual G4GeometryType  GetEntityType() const { return G4String("G4Trd"); }

                        // Visualisation functions

    void                DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent         GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4NURBS*      CreateNURBS        () const;

    void CheckAndSetAllParameters (G4double pdx1, G4double pdx2,
                             G4double pdy1, G4double pdy2,
                             G4double pdz);

    void SetAllParameters (G4double pdx1, G4double pdx2,
                             G4double pdy1, G4double pdy2,
                             G4double pdz);

protected:

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
    
    G4double fDx1,fDx2,fDy1,fDy2,fDz;

// Codes for faces (kPX=plus x face,kMY= minus y face etc)

  enum ESide {kUndefined, kPX,kMX,kPY,kMY,kPZ,kMZ};

};

#endif







