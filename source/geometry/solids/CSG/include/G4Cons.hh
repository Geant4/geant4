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
// $Id: G4Cons.hh,v 1.7 2002-10-28 11:43:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Cons
//
// Class description:
//
//   A G4Cons is, in the general case, a Phi segment of a cone, with
//   half-length fDz, inner and outer radii specified at -fDz and +fDz.
//   The Phi segment is described by a starting fSPhi angle, and the
//   +fDPhi delta angle for the shape.
//   If the delta angle is >=2*M_PI, the shape is treated as continuous
//   in Phi
//
//   Member Data:
//
//  fRmin1  inside radius at  -fDz
//  fRmin2  inside radius at  +fDz
//  fRmax1  outside radius at -fDz
//  fRmax2  outside radius at +fDz
//  fDz  half length in z
//
//  fSPhi  starting angle of the segment in radians
//  fDPhi  delta angle of the segment in radians
//
//   Note:
//      Internally fSPhi & fDPhi are adjusted so that fDPhi<=2PI,
//      and fDPhi+fSPhi<=2PI. This enables simpler comparisons to be
//      made with (say) Phi of a point.

// History:
// 19.3.94 P.Kent: Old C++ code converted to tolerant geometry
// 13.9.96 V.Grichine: Final modifications to commit
// --------------------------------------------------------------------

#ifndef G4Cons_HH
#define G4Cons_HH

#include "G4CSGSolid.hh"

class G4Cons : public G4CSGSolid
{
  public:  // with description

        G4Cons(const G4String& pName,
                     G4double pRmin1, G4double pRmax1,
                     G4double pRmin2, G4double pRmax2,
                     G4double pDz,
                     G4double pSPhi, G4double pDPhi);
         
        virtual ~G4Cons() ;

  // Accessors

        inline G4double    GetInnerRadiusMinusZ() const;
        inline G4double    GetOuterRadiusMinusZ() const;
        inline G4double    GetInnerRadiusPlusZ()  const;
        inline G4double    GetOuterRadiusPlusZ()  const;
  
        inline G4double    GetZHalfLength()       const;
  
        inline G4double    GetStartPhiAngle () const;
        inline G4double    GetDeltaPhiAngle () const;
  
  // Modifiers

        inline void    SetInnerRadiusMinusZ( G4double Rmin1 );
        inline void    SetOuterRadiusMinusZ( G4double Rmax1 );
        inline void    SetInnerRadiusPlusZ ( G4double Rmin2 );
        inline void    SetOuterRadiusPlusZ ( G4double Rmax2 );
         
        inline void    SetZHalfLength      ( G4double newDz );
        inline void    SetStartPhiAngle    ( G4double newSPhi);
        inline void    SetDeltaPhiAngle    ( G4double newDPhi);

  // Other methods for solid

        void ComputeDimensions(G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep);

        G4bool CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pmin, G4double& pmax) const;         

        EInside Inside(const G4ThreeVector& p) const;

        G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;

        G4double DistanceToIn(const G4ThreeVector& p,
                               const G4ThreeVector& v) const;
        G4double DistanceToIn(const G4ThreeVector& p) const;
        G4double DistanceToOut(const G4ThreeVector& p,
                               const G4ThreeVector& v,
                               const G4bool calcNorm=G4bool(false),
                                     G4bool *validNorm=0,
                                     G4ThreeVector *n=0) const;             
        G4double DistanceToOut(const G4ThreeVector& p) const;

        G4GeometryType  GetEntityType() const;

        G4std::ostream& StreamInfo(G4std::ostream& os) const;

  // Visualisation functions

        void          DescribeYourselfTo( G4VGraphicsScene& scene ) const;
        G4Polyhedron* CreatePolyhedron() const;
        G4NURBS*      CreateNURBS() const;

  public:  // without description
       
        //  Old access functions

        inline G4double    GetRmin1() const;
        inline G4double    GetRmax1() const;
        inline G4double    GetRmin2() const;
        inline G4double    GetRmax2() const;
  
        inline G4double    GetDz()    const;
  
        inline G4double    GetSPhi() const;
        inline G4double    GetDPhi() const;

  protected:
 
        G4ThreeVectorList*
        CreateRotatedVertices(const G4AffineTransform& pTransform) const;
  
        // Used by distanceToOut
  
        enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};
  
        // used by normal
  
        enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

  private:

        G4double fRmin1,fRmin2,
                 fRmax1,fRmax2,
                 fDz,
                 fSPhi,fDPhi;
};

#include "G4Cons.icc"

#endif
