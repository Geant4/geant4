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
// $Id: G4Sphere.hh,v 1.7 2002-10-28 11:43:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Sphere
//
// Class description:
//
//   A G4Sphere is, in the general case, section of a spherical shell,
//   between specified phii and theta angles
//
//   The phi and theta segments are described by a starting angle,
//   and the +ve delta angle for the shape.
//   If the delta angle is >=2*M_PI, or >=M_PI the shape is treated as
//   continuous in phi or theta respectively.
//
//   Theta must lie between 0-PI (incl).
//
//   Member Data:
//
//   fRmin  inner radius
//   fRmax  outer radius
//
//   fSPhi  starting angle of the segment in radians
//   fDPhi  delta angle of the segment in radians
//
//   fSTheta  starting angle of the segment in radians
//   fDTheta  delta angle of the segment in radians
//
//     
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.

// History:
// 28.3.94 P.Kent: old C++ code converted to tolerant geometry
// 17.9.96 V.Grichine: final modifications to commit
// --------------------------------------------------------------------

#ifndef G4Sphere_HH
#define G4Sphere_HH

#include "G4CSGSolid.hh"

class G4Sphere : public G4CSGSolid
{
  public:  // with description

    G4Sphere(const G4String& pName,
                   G4double pRmin, G4double pRmax,
                   G4double pSPhi, G4double pDPhi,
                   G4double pSTheta, G4double pDTheta);
       
    virtual ~G4Sphere() ;
    
    // Accessors
       
    inline G4double GetInsideRadius   () const;
    inline G4double GetOuterRadius    () const;
    inline G4double GetStartPhiAngle  () const;
    inline G4double GetDeltaPhiAngle  () const;
    inline G4double GetStartThetaAngle() const;
    inline G4double GetDeltaThetaAngle() const;

    // Modifiers

    inline void SetInsideRadius   (G4double newRmin);
    inline void SetOuterRadius    (G4double newRmax);
    inline void SetStartPhiAngle  (G4double newSphi);
    inline void SetDeltaPhiAngle  (G4double newDphi);
    inline void SetStartThetaAngle(G4double newSTheta);
    inline void SetDeltaThetaAngle(G4double newDTheta);

    // Methods for solid

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

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
 
    G4std::ostream& StreamInfo(G4std::ostream& os) const;

    // Visualisation functions
  
    void          DescribeYourselfTo(G4VGraphicsScene& scene) const;
    G4Polyhedron* CreatePolyhedron() const;
    G4NURBS*      CreateNURBS() const;

  public:  // without description
   
    // Old access functions

    inline G4double  GetRmin()   const;
    inline G4double  GetRmax()   const;
    inline G4double  GetSPhi()   const;
    inline G4double  GetDPhi()   const;
    inline G4double  GetSTheta() const;
    inline G4double  GetDTheta() const;

  protected:
 
    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform,
                                G4int& noPolygonVertices) const;
  
    // Used by distanceToOut
  
    enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kSTheta,kETheta};
  
    // used by normal
  
    enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNSTheta,kNETheta};

  private:

    G4double fRmin,fRmax,
             fSPhi,fDPhi,
             fSTheta,fDTheta;
};

#include "G4Sphere.icc"

#endif
