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
// $Id: G4Box.hh,v 1.8 2002-10-28 11:43:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4Box
//
// Class description:
//
//   A Box is a cuboid of given half lengths dx,dy,dz. The Box is
//   centred on the origin with sides parallel to the x/y/z axes.

// History:
// 30.06.95 P.Kent: Converted from source code developed end 94
// 27.03.96 J.Allison: Added virtual functions DescribeYourselfTo() and
//                     SendWireframeTo(G4VGraphicsModel&)
// 22.07.96 J.Allison: Changed G4VGraphicsModel to G4VGraphicsScene
//                     and SendPolyhedronTo() to CreatePolyhedron()
// 27.03.98 J.Apostolakis: Inherit from G4CSGSolid (not G4VSolid)
// 18.11.99 J.Apostolakis, V.Grichine: kUndefined was added to ESide
// --------------------------------------------------------------------

#ifndef G4BOX_HH
#define G4BOX_HH

#include "G4CSGSolid.hh"

class G4Box : public G4CSGSolid 
{
  public:  // with description

    G4Box(const G4String& pName, G4double pX, G4double pY, G4double pZ);
      // Construct a box with name, and half lengths pX,pY,pZ

    virtual ~G4Box();

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4bool CalculateExtent(const EAxis pAxis,
			   const G4VoxelLimits& pVoxelLimit,
			   const G4AffineTransform& pTransform,
			         G4double& pmin, G4double& pmax) const;

  // Accessors and modifiers

    inline G4double GetXHalfLength() const;
    
    inline G4double GetYHalfLength() const;
    
    inline G4double GetZHalfLength() const;

    void SetXHalfLength(G4double dx) ;

    void SetYHalfLength(G4double dy) ;

    void SetZHalfLength(G4double dz) ;

  // Methods for solid

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;

    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
			   const G4bool calcNorm=false,
			         G4bool *validNorm=0, G4ThreeVector *n=0) const;

    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4std::ostream& StreamInfo(G4std::ostream& os) const;

  // Functions for visualization

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4NURBS*      CreateNURBS        () const;

  protected:  // with description

    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

  protected:  // without description

    enum ESide {kUndefined,kPX,kMX,kPY,kMY,kPZ,kMZ};
      // Codes for faces (kPX=plus x face,kMY= minus y face etc)

  private:

    G4double fDx,fDy,fDz;
};

#include "G4Box.icc"

#endif
