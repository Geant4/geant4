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
// $Id: G4TwistedGenTrap.hh,v 1.2 2005-03-07 13:24:06 link Exp $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedGenTrap
//
// Class description:
//
//  G4TwistedTrap is a general trapezoid defined with caps of different shape
//  and size, an twisted along the Z axis.

// Author:
//                 O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#ifndef __G4TWISTEDGENTRAP__
#define __G4TWISTEDGENTRAP__

#include "G4VSolid.hh"
#include "G4TwistedTrapAlphaSide.hh"
#include "G4TwistedTrapParallelSide.hh"
#include "G4FlatTrapSide.hh" 

class G4SolidExtentList;
class G4ClippablePolygon;

class G4TwistedGenTrap: public G4VSolid
{
 public:  // with description
 
  G4TwistedGenTrap(const G4String &pname,         // Name of instance
		   G4double      PhiTwist,    // twist angle
		   G4double      pDz,         // half z lenght
		   G4double      pTheta,      // direction between end planes
		   G4double      pPhi,        // defined by polar and azimutal angles.
		   G4double      pDy1,        // half y length at -pDz
		   G4double      pDx1,        // half x length at -pDz,-pDy
		   G4double      pDx2,        // half x length at -pDz,+pDy
		   G4double      pDy2,        // half y length at +pDz
		   G4double      pDx3,        // half x length at +pDz,-pDy
		   G4double      pDx4,        // half x length at +pDz,+pDy
		   G4double      pAlph        // tilt angle at +pDz
		   );
  
  virtual ~G4TwistedGenTrap();
             
  void ComputeDimensions(G4VPVParameterisation*    ,
                         const G4int               ,
                         const G4VPhysicalVolume*  );
 
  G4bool CalculateExtent(const EAxis               pAxis,
                         const G4VoxelLimits      &pVoxelLimit,
                         const G4AffineTransform  &pTransform,
                               G4double           &pMin,
                               G4double           &pMax ) const;

  G4double DistanceToIn (const G4ThreeVector &p,
                         const G4ThreeVector &v ) const;

  G4double DistanceToIn (const G4ThreeVector &p ) const;
   
  G4double DistanceToOut(const G4ThreeVector &p, 
                         const G4ThreeVector &v,
                         const G4bool         calcnorm=G4bool(false),
                               G4bool        *validnorm=0, 
                               G4ThreeVector *n=0 ) const;

  G4double DistanceToOut(const G4ThreeVector &p) const;
  
  EInside Inside (const G4ThreeVector &p) const;

  G4ThreeVector SurfaceNormal(const G4ThreeVector &p) const;

  inline G4double GetCubicVolume() ;

  void            DescribeYourselfTo (G4VGraphicsScene &scene) const;
  G4Polyhedron   *CreatePolyhedron   () const;
  G4NURBS        *CreateNURBS        () const;

  std::ostream &StreamInfo(std::ostream& os) const;

  // accessors
  
  inline G4double GetPhiTwist    () const { return fPhiTwist   ; }

  inline G4double GetDx1   () const { return fDx1 ; } 
  inline G4double GetDx2   () const { return fDx2 ; } 
  inline G4double GetDx3   () const { return fDx3 ; } 
  inline G4double GetDx4   () const { return fDx4 ; } 
  inline G4double GetDy1    () const { return fDy1 ; } 
  inline G4double GetDy2    () const { return fDy2 ; } 
  inline G4double GetDz () const { return fDz; }
  inline G4double GetPhi () const { return fPhi ; }
  inline G4double GetTheta () const { return fTheta ; }
  inline G4double GetAlpha () const { return fAlph ; }

  inline G4double xAxisMin(G4double phi) const ;
  inline G4double xAxisMax(G4double phi) const ;
  inline G4double yAxisMax(G4double phi) const ;

  G4VisExtent     GetExtent    () const;
  G4GeometryType  GetEntityType() const;

 protected:  // with description

  G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.


 public:  // without description

 private:
 
  void         CreateSurfaces();

 private:
 
  G4double fDx ;  // temp

  G4double fTheta;   
  G4double fPhi ;

  G4double fDy1;   
  G4double fDx1;     
  G4double fDx2;     
  
  G4double fDy2;   
  G4double fDx3;     
  G4double fDx4;     

  G4double fDz;         // Half-length along the z axis
  
  G4double fAlph ;
  G4double fTAlph ;    // tan(fAlph)

  G4double fdeltaX ;
  G4double fdeltaY ;
    
  G4double fPhiTwist;   // twist angle ( dphi in surface equation)

  G4double fAngleSide;
     
  G4VSurface *fLowerEndcap ;  // surface of -ve z
  G4VSurface *fUpperEndcap ;  // surface of +ve z
  
  G4VSurface *fSide0 ;         // Twisted Side at phi = 0 deg
  G4VSurface *fSide90 ;        // Twisted Side at phi = 90 deg
  G4VSurface *fSide180 ;       // Twisted Side at phi = 180 deg
  G4VSurface *fSide270 ;       // Twisted Side at phi = 270 deg

  G4double fCubicVolume ;      // volume of the twisted trapezoid

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
G4double G4TwistedGenTrap::GetCubicVolume()
{

  
  if(fCubicVolume != 0.) ;
  else   fCubicVolume = 1 ;
  return fCubicVolume;
}


inline 
G4double G4TwistedGenTrap::xAxisMin(G4double phi) const
{
  return ( - (2*fDx2*(fPhiTwist - 2*phi) + 2*fDx4*(fPhiTwist + 2*phi) + (2*fDy1*(fPhiTwist - 2*phi) + 2*fDy2*(fPhiTwist + 2*phi))*fTAlph)/(4.*fPhiTwist) ) ;
}

inline
G4double G4TwistedGenTrap::xAxisMax(G4double phi) const
{
  return   (2*fDx2*(fPhiTwist - 2*phi) + 2*fDx4*(fPhiTwist + 2*phi) - (2*fDy1*(fPhiTwist - 2*phi) + 2*fDy2*(fPhiTwist + 2*phi))*fTAlph)/(4.*fPhiTwist) ;
}

inline 
G4double G4TwistedGenTrap::yAxisMax(G4double phi) const
{
  return ( fDy1 + fDy2 + ( fDy2 - fDy1  )  * ( 2 * phi ) / fPhiTwist ) ;
}






#endif
