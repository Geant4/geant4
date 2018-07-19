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
// $Id: G4TwistedTubs.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTubs
//
// Class description:
//
//  G4TwistedTubs is a sort of twisted cylinder.
//  A twisted cylinder which is placed along with z-axis and is
//  separated into phi-segments should become a hyperboloid, and
//  its each segmented piece should be tilted with a stereo angle. 
//  G4TwistedTubs is a G4VSolid.
//  It can have inner & outer surfaces as well as G4TwistedTubs, 
//  but cannot has different stereo angles between the inner surface
//  and outer surface.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTUBS__
#define __G4TWISTEDTUBS__

#include "G4VSolid.hh"
#include "G4TwistTubsFlatSide.hh"
#include "G4TwistTubsSide.hh"
#include "G4TwistTubsHypeSide.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

class G4TwistedTubs : public G4VSolid
{
 public:  // with description
 
  G4TwistedTubs(const G4String &pname,         // Name of instance
                      G4double  twistedangle,  // Twisted angle
                      G4double  endinnerrad,   // Inner radius at endcap 
                      G4double  endouterrad,   // Outer radius at endcap 
                      G4double  halfzlen,      // half z length 
                      G4double  dphi);         // Phi angle of a segment
                      
  G4TwistedTubs(const G4String &pname,         // Name of instance
                      G4double  twistedangle,  // Stereo angle
                      G4double  endinnerrad,   // Inner radius at endcap 
                      G4double  endouterrad,   // Outer radius at endcap 
                      G4double  halfzlen,      // half z length 
                      G4int     nseg,          // Number of segments in totalPhi
                      G4double  totphi);       // Total angle of all segments
                      
  G4TwistedTubs(const G4String &pname,         // Name of instance
                      G4double  twistedangle,  // Twisted angle
                      G4double  innerrad,      // Inner radius at z=0 
                      G4double  outerrad,      // Outer radius at z=0 
                      G4double  negativeEndz,  // -ve z endplate
                      G4double  positiveEndz,  // +ve z endplate
                      G4double  dphi);         // Phi angle of a segment

  G4TwistedTubs(const G4String &pname,         // Name of instance
                      G4double  twistedangle,  // Stereo angle
                      G4double  innerrad,      // Inner radius at z=0 
                      G4double  outerrad,      // Outer radius at z=0 
                      G4double  negativeEndz,  // -ve z endplate
                      G4double  positiveEndz,  // +ve z endplate
                      G4int     nseg,          // Number of segments in totalPhi
                      G4double  totphi);       // Total angle of all segments

  virtual ~G4TwistedTubs();
             
  void ComputeDimensions(G4VPVParameterisation   *  /* p  */ ,
                         const G4int                /* n  */ ,
                         const G4VPhysicalVolume *  /* prep */ );

  void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const; 

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

  void            DescribeYourselfTo (G4VGraphicsScene &scene) const;
  G4Polyhedron   *CreatePolyhedron   () const;
  G4Polyhedron   *GetPolyhedron      () const;

  std::ostream &StreamInfo(std::ostream& os) const;

  // accessors
  
  inline G4double GetDPhi        () const { return fDPhi       ; }
  inline G4double GetPhiTwist    () const { return fPhiTwist   ; }
  inline G4double GetInnerRadius () const { return fInnerRadius; }
  inline G4double GetOuterRadius () const { return fOuterRadius; }
  inline G4double GetInnerStereo () const { return fInnerStereo; }
  inline G4double GetOuterStereo () const { return fOuterStereo; }
  inline G4double GetZHalfLength () const { return fZHalfLength; }
  inline G4double GetKappa       () const { return fKappa      ; }

  inline G4double GetTanInnerStereo () const { return fTanInnerStereo  ; }
  inline G4double GetTanInnerStereo2() const { return fTanInnerStereo2 ; }
  inline G4double GetTanOuterStereo () const { return fTanOuterStereo  ; }
  inline G4double GetTanOuterStereo2() const { return fTanOuterStereo2 ; }
  
  inline G4double GetEndZ           (G4int i) const { return fEndZ[i]  ; }
  inline G4double GetEndPhi         (G4int i) const { return fEndPhi[i]; }
  inline G4double GetEndInnerRadius (G4int i) const
                  { return fEndInnerRadius[i]; }
  inline G4double GetEndOuterRadius (G4int i) const
                  { return fEndOuterRadius[i]; }
  inline G4double GetEndInnerRadius () const 
                  { return (fEndInnerRadius[0] > fEndInnerRadius[1] ?
                    fEndInnerRadius[0] : fEndInnerRadius[1]); }
  inline G4double GetEndOuterRadius () const
                  { return (fEndOuterRadius[0] > fEndOuterRadius[1] ?
                    fEndOuterRadius[0] : fEndOuterRadius[1]); }
  
  G4VisExtent     GetExtent    () const;
  G4GeometryType  GetEntityType() const;
  G4VSolid* Clone() const;

  G4double GetCubicVolume();
    // Returns an estimation of the geometrical cubic volume of the
    // solid. Caches the computed value once computed the first time.
  G4double GetSurfaceArea();
    // Returns an estimation of the geometrical surface area of the
    // solid. Caches the computed value once computed the first time.

  G4ThreeVector GetPointOnSurface() const ;

 public:  // without description

  G4TwistedTubs(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4TwistedTubs(const G4TwistedTubs& rhs);
  G4TwistedTubs& operator=(const G4TwistedTubs& rhs); 
    // Copy constructor and assignment operator.

#ifdef G4TWISTDEBUG
  G4VTwistSurface * GetOuterHype() const { return fOuterHype; }
#endif
  
 private:
 
  inline void  SetFields(G4double phitwist, G4double innerrad,
                         G4double outerrad,
                         G4double negativeEndz, G4double positiveEndz);
                     
  void         CreateSurfaces();

 private:
 
  G4double fPhiTwist;       // Twist angle from -fZHalfLength to fZHalfLength
  G4double fInnerRadius;    // Inner-hype radius at z=0
  G4double fOuterRadius;    // Outer-hype radius at z=0
  G4double fEndZ[2];        // z at endcaps, [0] = -ve z, [1] = +ve z
  G4double fDPhi;           // Phi-width of a segment fDPhi > 0
  G4double fZHalfLength;    // Half length along z-axis
 
  G4double fInnerStereo;       // Inner-hype stereo angle
  G4double fOuterStereo;       // Outer-hype stereo angle
  G4double fTanInnerStereo;    // std::tan(innerStereoAngle)
  G4double fTanOuterStereo;    // std::tan(outerStereoAngle)
  G4double fKappa;             // std::tan(fPhiTwist/2)/fZHalfLen;
  G4double fEndInnerRadius[2]; // Inner-hype radii endcaps [0] -ve z, [1] +ve z
  G4double fEndOuterRadius[2]; // Outer-hype radii endcaps [0] -ve z, [1] +ve z
  G4double fEndPhi[2];         // Phi endcaps, [0] = -ve z, [1] = +ve z
  
  G4double fInnerRadius2;      // fInnerRadius * fInnerRadius
  G4double fOuterRadius2;      // fOuterRadius * fOuterRadius
  G4double fTanInnerStereo2;   // fInnerRadius * fInnerRadius
  G4double fTanOuterStereo2;   // fInnerRadius * fInnerRadius
  G4double fEndZ2[2];          // fEndZ * fEndZ
  
  G4VTwistSurface *fLowerEndcap;    // Surface of -ve z
  G4VTwistSurface *fUpperEndcap;    // Surface of +ve z
  G4VTwistSurface *fLatterTwisted;  // Surface of -ve phi
  G4VTwistSurface *fFormerTwisted;  // Surface of +ve phi
  G4VTwistSurface *fInnerHype;      // Surface of -ve r
  G4VTwistSurface *fOuterHype;      // Surface of +ve r

  G4double fCubicVolume;       // Cached value for cubic volume
  G4double fSurfaceArea;       // Cached value for surface area

  mutable G4bool fRebuildPolyhedron;
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

//---------------------
// inline functions
//---------------------

inline
void G4TwistedTubs::SetFields(G4double phitwist, G4double innerrad, 
                              G4double outerrad, G4double negativeEndz, 
                              G4double positiveEndz)
{
   fCubicVolume  = 0.;
   fPhiTwist     = phitwist;
   fEndZ[0]      = negativeEndz;
   fEndZ[1]      = positiveEndz;
   fEndZ2[0]     = fEndZ[0] * fEndZ[0];
   fEndZ2[1]     = fEndZ[1] * fEndZ[1];
   fInnerRadius  = innerrad;
   fOuterRadius  = outerrad;
   fInnerRadius2 = fInnerRadius * fInnerRadius;
   fOuterRadius2 = fOuterRadius * fOuterRadius;
   
   if (std::fabs(fEndZ[0]) >= std::fabs(fEndZ[1])) {
      fZHalfLength = std::fabs(fEndZ[0]);
   } else {
      fZHalfLength = std::fabs(fEndZ[1]);
   }

   G4double parity         = (fPhiTwist > 0 ? 1 : -1); 
   G4double tanHalfTwist   = std::tan(0.5 * fPhiTwist);
   G4double innerNumerator = std::fabs(fInnerRadius * tanHalfTwist) * parity;
   G4double outerNumerator = std::fabs(fOuterRadius * tanHalfTwist) * parity;

   fTanInnerStereo    = innerNumerator / fZHalfLength; 
   fTanOuterStereo    = outerNumerator / fZHalfLength; 
   fTanInnerStereo2   = fTanInnerStereo * fTanInnerStereo;
   fTanOuterStereo2   = fTanOuterStereo * fTanOuterStereo;
   fInnerStereo       = std::atan2(innerNumerator,  fZHalfLength); 
   fOuterStereo       = std::atan2(outerNumerator,  fZHalfLength); 
   fEndInnerRadius[0] = std::sqrt(fInnerRadius2 + fEndZ2[0] * fTanInnerStereo2);
   fEndInnerRadius[1] = std::sqrt(fInnerRadius2 + fEndZ2[1] * fTanInnerStereo2);
   fEndOuterRadius[0] = std::sqrt(fOuterRadius2 + fEndZ2[0] * fTanOuterStereo2);
   fEndOuterRadius[1] = std::sqrt(fOuterRadius2 + fEndZ2[1] * fTanOuterStereo2);

   fKappa          = tanHalfTwist / fZHalfLength;
   fEndPhi[0]      = std::atan2(fEndZ[0] * tanHalfTwist, fZHalfLength);
   fEndPhi[1]      = std::atan2(fEndZ[1] * tanHalfTwist, fZHalfLength);

#ifdef G4TWISTDEBUG
   G4cout << "/********* G4TwistedTubs::SetFields() Field Parameters ***************** " << G4endl;
   G4cout << "/*   fPhiTwist                  : " << fPhiTwist << G4endl;
   G4cout << "/*   fEndZ(0, 1)                : " << fEndZ[0] << " , " << fEndZ[1] << G4endl; 
   G4cout << "/*   fEndPhi(0, 1)              : " << fEndPhi[0] << " , " << fEndPhi[1] << G4endl; 
   G4cout << "/*   fInnerRadius, fOuterRadius : " << fInnerRadius << " , " << fOuterRadius << G4endl; 
   G4cout << "/*   fEndInnerRadius(0, 1)      : " << fEndInnerRadius[0] << " , " 
          << fEndInnerRadius[1] << G4endl; 
   G4cout << "/*   fEndOuterRadius(0, 1)      : " << fEndOuterRadius[0] << " , " 
          << fEndOuterRadius[1] << G4endl; 
   G4cout << "/*   fInnerStereo, fOuterStereo : " << fInnerStereo << " , " << fOuterStereo << G4endl; 
   G4cout << "/*   tanHalfTwist, fKappa       : " << tanHalfTwist << " , " << fKappa << G4endl; 
   G4cout << "/*********************************************************************** " << G4endl;
#endif
}

#endif
