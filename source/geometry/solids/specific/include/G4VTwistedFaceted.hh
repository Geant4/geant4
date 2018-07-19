//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VTwistedFaceted.hh 104316 2017-05-24 13:04:23Z gcosmo $
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

// Author:
//
//   27-Oct-2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#ifndef __G4VTWISTEDFACETED__
#define __G4VTWISTEDFACETED__

#include "G4VSolid.hh"
#include "G4TwistTrapAlphaSide.hh"
#include "G4TwistTrapParallelSide.hh"
#include "G4TwistBoxSide.hh"
#include "G4TwistTrapFlatSide.hh" 

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
 
  virtual void BoundingLimits(G4ThreeVector &pMin, G4ThreeVector &pMax) const;

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

  G4ThreeVector GetPointOnSurface() const;
  G4ThreeVector GetPointInSolid(G4double z) const;
  
  virtual inline G4double GetCubicVolume() ;
  virtual inline G4double GetSurfaceArea() ;

  virtual void            DescribeYourselfTo (G4VGraphicsScene &scene) const;
  virtual G4Polyhedron   *CreatePolyhedron   () const ;
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

 public:  // without description

  G4VTwistedFaceted(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4VTwistedFaceted(const G4VTwistedFaceted& rhs);
  G4VTwistedFaceted& operator=(const G4VTwistedFaceted& rhs); 
    // Copy constructor and assignment operator.

 protected:  // with description

  mutable G4bool fRebuildPolyhedron;
  mutable G4Polyhedron* fpPolyhedron;  // pointer to polyhedron for vis

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

  G4VTwistSurface *fLowerEndcap ;  // surface of -ve z
  G4VTwistSurface *fUpperEndcap ;  // surface of +ve z
  
  G4VTwistSurface *fSide0 ;         // Twisted Side at phi = 0 deg
  G4VTwistSurface *fSide90 ;        // Twisted Side at phi = 90 deg
  G4VTwistSurface *fSide180 ;       // Twisted Side at phi = 180 deg
  G4VTwistSurface *fSide270 ;       // Twisted Side at phi = 270 deg

  G4double fCubicVolume ;      // volume of the solid
  G4double fSurfaceArea ;      // area of the solid

  class LastState              // last Inside result
  {
    public:
      LastState()
      {
        p.set(kInfinity,kInfinity,kInfinity); inside = kOutside;
      }
      ~LastState(){}
      LastState(const LastState& r) : p(r.p), inside(r.inside){}
      LastState& operator=(const LastState& r)
      {
        if (this == &r)  { return *this; }
        p = r.p; inside = r.inside;
        return *this;
      }
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
        surface = new G4VTwistSurface*[1];
      }
      ~LastVector()
      {
        delete [] surface;
      }
      LastVector(const LastVector& r) : p(r.p), vec(r.vec)
      {
        surface = new G4VTwistSurface*[1];
        surface[0] = r.surface[0];
      }
      LastVector& operator=(const LastVector& r)
      {
        if (&r == this)  { return *this; }
        p = r.p; vec = r.vec;
        delete [] surface; surface = new G4VTwistSurface*[1];
        surface[0] = r.surface[0];
        return *this;
      }
    public:
      G4ThreeVector   p;
      G4ThreeVector   vec;
      G4VTwistSurface **surface;
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
      LastValue(const LastValue& r) : p(r.p), value(r.value){}
      LastValue& operator=(const LastValue& r)
      {
        if (this == &r)  { return *this; }
        p = r.p; value = r.value;
        return *this;
      }
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
      LastValueWithDoubleVector(const LastValueWithDoubleVector& r)
        : p(r.p), vec(r.vec), value(r.value){}
      LastValueWithDoubleVector& operator=(const LastValueWithDoubleVector& r)
      {
        if (this == &r)  { return *this; }
        p = r.p; vec = r.vec; value = r.value;
        return *this;
      }
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
G4double G4VTwistedFaceted::GetSurfaceArea()
{
  if(fSurfaceArea != 0.) ;
  else   fSurfaceArea = G4VSolid::GetSurfaceArea();
  return fSurfaceArea;
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
    - u*( ( GetValueD(phi)-GetValueA(phi) ) / ( 2 * GetValueB(phi) ) - ftg );
}

#endif
