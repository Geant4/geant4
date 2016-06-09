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
// $Id: G4VTwistedFaceted.hh,v 1.3 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4VTwistedFaceted
//
// Class description:
//
//  G4VTwistedFaceted is an abstract base class for twisted boxoids:
//  G4TwistedTrd, G4TwistedTrap and G4TwistedBox

// Author:  O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#ifndef __G4VTWISTEDFACETED__
#define __G4VTWISTEDFACETED__

#include "G4VSolid.hh"
#include "G4TwistedTrapAlphaSide.hh"
#include "G4TwistedTrapParallelSide.hh"
#include "G4TwistedTrapBoxSide.hh"
#include "G4FlatTrapSide.hh" 

class G4SolidExtentList;
class G4ClippablePolygon;

class G4VTwistedFaceted: public G4VSolid
{
 public:  // with description
 
  G4VTwistedFaceted(const G4String &pname,    // Name of instance
                          G4double PhiTwist,  // twist angle
                          G4double pDz,       // half z lenght
                          G4double pTheta,  // direction between end planes
                          G4double pPhi,    // defined by polar & azim. angles
                          G4double pDy1,    // half y length at -pDz
                          G4double pDx1,    // half x length at -pDz,-pDy
                          G4double pDx2,    // half x length at -pDz,+pDy
                          G4double pDy2,    // half y length at +pDz
                          G4double pDx3,    // half x length at +pDz,-pDy
                          G4double pDx4,    // half x length at +pDz,+pDy
                          G4double pAlph    // tilt angle at +pDz
                   );
  
  virtual ~G4VTwistedFaceted();
             
  virtual void ComputeDimensions(G4VPVParameterisation*,
                                 const G4int,
                                 const G4VPhysicalVolume*  );
 
  virtual G4bool CalculateExtent(const EAxis               pAxis,
                                 const G4VoxelLimits      &pVoxelLimit,
                                 const G4AffineTransform  &pTransform,
                                       G4double           &pMin,
                                       G4double           &pMax ) const;

  virtual G4double DistanceToIn (const G4ThreeVector &p,
                                 const G4ThreeVector &v ) const;

  virtual G4double DistanceToIn (const G4ThreeVector &p ) const;
   
  virtual G4double DistanceToOut(const G4ThreeVector &p, 
                                 const G4ThreeVector &v,
                                 const G4bool         calcnorm  = false,
                                       G4bool        *validnorm = 0, 
                                       G4ThreeVector *n=0 ) const;

  virtual G4double DistanceToOut(const G4ThreeVector &p) const;
  
  virtual EInside Inside (const G4ThreeVector &p) const;

  virtual G4ThreeVector SurfaceNormal(const G4ThreeVector &p) const;

  virtual inline G4double GetCubicVolume() ;

  virtual void            DescribeYourselfTo (G4VGraphicsScene &scene) const;
  virtual G4Polyhedron   *CreatePolyhedron   () const = 0 ;
  virtual G4NURBS        *CreateNURBS        () const;
  virtual G4Polyhedron   *GetPolyhedron      () const;

  virtual std::ostream &StreamInfo(std::ostream& os) const;

  // accessors
  
  inline G4double GetTwistAngle    () const { return fPhiTwist; }

  inline G4double GetDx1   () const { return fDx1   ; } 
  inline G4double GetDx2   () const { return fDx2   ; } 
  inline G4double GetDx3   () const { return fDx3   ; } 
  inline G4double GetDx4   () const { return fDx4   ; } 
  inline G4double GetDy1   () const { return fDy1   ; } 
  inline G4double GetDy2   () const { return fDy2   ; } 
  inline G4double GetDz    () const { return fDz    ; }
  inline G4double GetPhi   () const { return fPhi   ; }
  inline G4double GetTheta () const { return fTheta ; }
  inline G4double GetAlpha () const { return fAlph  ; }

  inline G4double Xcoef(G4double u,G4double phi, G4double ftg) const ;
    // For calculating the w(u) function

  inline G4double GetValueA(G4double phi) const;
  inline G4double GetValueB(G4double phi) const;
  inline G4double GetValueD(G4double phi) const;

  virtual G4VisExtent     GetExtent    () const;
  virtual G4GeometryType  GetEntityType() const;

 protected:  // with description

  G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

 private:
 
  void CreateSurfaces();

 private:
 
  G4double fTheta;   
  G4double fPhi ;

  G4double fDy1;   
  G4double fDx1;     
  G4double fDx2;     
  
  G4double fDy2;   
  G4double fDx3;     
  G4double fDx4;     

  G4double fDz;        // Half-length along the z axis

  G4double fDx ;       // maximum side in x 
  G4double fDy ;       // maximum side in y

  G4double fAlph ;
  G4double fTAlph ;    // std::tan(fAlph)

  G4double fdeltaX ;
  G4double fdeltaY ;
    
  G4double fPhiTwist;  // twist angle ( dphi in surface equation)

  G4double fAngleSide;
     
  G4VSurface *fLowerEndcap ;  // surface of -ve z
  G4VSurface *fUpperEndcap ;  // surface of +ve z
  
  G4VSurface *fSide0 ;         // Twisted Side at phi = 0 deg
  G4VSurface *fSide90 ;        // Twisted Side at phi = 90 deg
  G4VSurface *fSide180 ;       // Twisted Side at phi = 180 deg
  G4VSurface *fSide270 ;       // Twisted Side at phi = 270 deg

  G4double fCubicVolume ;      // volume of the twisted trapezoid

  mutable G4Polyhedron* fpPolyhedron;  // pointer to polyhedron for vis

  class LastState              // last Inside result
  {
    public:
      LastState()
      {
        p.set(kInfinity,kInfinity,kInfinity);
        inside = kOutside;
      }
      ~LastState(){}
    public:
      G4ThreeVector p;
        EInside       inside;
  };
              
  class LastVector             // last SurfaceNormal result
  {
    public:
      LastVector()
      {
        p.set(kInfinity,kInfinity,kInfinity);
        vec.set(kInfinity,kInfinity,kInfinity);
        surface = new G4VSurface*[1];
      }
      ~LastVector()
      {
        delete [] surface;
      }
    public:
      G4ThreeVector   p;
      G4ThreeVector   vec;
      G4VSurface    **surface;
  };

  class LastValue              // last G4double value
  {
    public:
      LastValue()
      {
        p.set(kInfinity,kInfinity,kInfinity);
        value = DBL_MAX;
      }
      ~LastValue(){}
    public:
      G4ThreeVector p;
      G4double      value;
  };
              
  class LastValueWithDoubleVector   // last G4double value
  {
    public:
      LastValueWithDoubleVector()
      {
        p.set(kInfinity,kInfinity,kInfinity);
        vec.set(kInfinity,kInfinity,kInfinity);
        value = DBL_MAX;
      }
      ~LastValueWithDoubleVector(){}
      public:
        G4ThreeVector p;
        G4ThreeVector vec;
        G4double      value;
  };
              
  LastState    fLastInside;
  LastVector   fLastNormal;
  LastValue    fLastDistanceToIn;
  LastValue    fLastDistanceToOut;
  LastValueWithDoubleVector   fLastDistanceToInWithV;
  LastValueWithDoubleVector   fLastDistanceToOutWithV;

 };

//=====================================================================

inline
G4double G4VTwistedFaceted::GetCubicVolume()
{

  
  if(fCubicVolume != 0.) ;
  else   fCubicVolume = 2 * fDz
                      * ( ( fDx1 + fDx2 ) * fDy1 + ( fDx3 + fDx4 ) * fDy2  );
  return fCubicVolume;
}

inline
G4double G4VTwistedFaceted::GetValueA(G4double phi) const
{
  return ( fDx4 + fDx2  + ( fDx4 - fDx2 ) * ( 2 * phi ) / fPhiTwist  ) ;
}

inline
G4double G4VTwistedFaceted::GetValueD(G4double phi) const
{
  return ( fDx3 + fDx1  + ( fDx3 - fDx1 ) * ( 2 * phi ) / fPhiTwist  ) ;
} 

inline 
G4double G4VTwistedFaceted::GetValueB(G4double phi) const
{
  return ( fDy2 + fDy1  + ( fDy2 - fDy1 ) * ( 2 * phi ) / fPhiTwist ) ;
}

inline
G4double G4VTwistedFaceted::Xcoef(G4double u, G4double phi, G4double ftg) const 
{
  return GetValueA(phi)/2. + (GetValueD(phi)-GetValueA(phi))/4. 
    - u*( ( GetValueD(phi)-GetValueA(phi) ) / ( 2 * GetValueB(phi) ) + ftg );
}

#endif
