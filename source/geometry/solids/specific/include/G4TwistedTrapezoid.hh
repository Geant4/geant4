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
// $Id: G4TwistedTrapezoid.hh,v 1.2 2004-08-27 13:37:12 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrapezoid
//
// Class description:
//
//  G4TwistedTrapezoid is a sort of twisted box.
//
// Author:
//                 O.Link (Oliver.Link@cern.ch)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
//                 see: Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
// --------------------------------------------------------------------

#ifndef __G4TWISTEDTRAPEZOID__
#define __G4TWISTEDTRAPEZOID__

#include "G4VSolid.hh"
#include "G4TwistedTrapSide.hh"

class G4SolidExtentList;
class G4ClippablePolygon;

class G4TwistedTrapezoid : public G4VSolid
{
 public:  // with description
 
  G4TwistedTrapezoid(const G4String &pname,         // Name of instance
		     G4double  twistedangle,  // Twisted angle
		     G4double  halfzlen,      // half z length 
		     G4double  halfSideX,         // half length in x
		     G4double  halfSideY,         // half length in y
                      G4double  dphi);         // Phi angle of a segment
                      
   virtual ~G4TwistedTrapezoid();
             
  void ComputeDimensions(G4VPVParameterisation   *  /* p  */ ,
                         const G4int                /* n  */ ,
                         const G4VPhysicalVolume *  /* prep */ );
 
  G4bool CalculateExtent(const EAxis               paxis,
                         const G4VoxelLimits      &pvoxellimit,
                         const G4AffineTransform  &ptransform,
                               G4double           &pmin,
                               G4double           &pmax ) const;

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
  G4NURBS        *CreateNURBS        () const;

  std::ostream &StreamInfo(std::ostream& os) const;

  // accessors
  
  inline G4double GetPhiTwist    () const { return fPhiTwist   ; }
  inline G4double GetZHalfLength () const { return fZHalfLength; }

  inline G4double GetHalfSideX   () const { return fHalfSides[0] ; } 
  inline G4double GetHalfSideY   () const { return fHalfSides[1] ; } 
  
  G4VisExtent     GetExtent    () const;
  G4GeometryType  GetEntityType() const;

 public:  // without description

 private:
 
  inline void  SetFields(G4double phitwist, G4double fHalfZ, 
			 G4double fHalfSideX, G4double fHalfSideY);
                     
  void         CreateSurfaces();

  
  static void  AddPolyToExtent( const G4ThreeVector     &v0,
                                const G4ThreeVector     &v1,
                                const G4ThreeVector     &w1,
                                const G4ThreeVector     &w0,
                                const G4VoxelLimits     &voxellimit,
                                const EAxis              axis,
                                      G4SolidExtentList &extentlist );
 private:
 
  G4double fPhiTwist;       // Twist angle from -fZHalfLength to fZHalfLength
  G4double fZHalfLength;    // Half length along z-axis

  G4double fHalfSides[2];   // Half length along x and y axis
     
  G4VSurface *fLowerEndcap ;  // surface of -ve z
  G4VSurface *fUpperEndcap ;  // surface of +ve z
  
  G4VSurface *fSide0 ;         // Twisted Side at phi = 0 deg
  G4VSurface *fSide90 ;        // Twisted Side at phi = 90 deg
  G4VSurface *fSide180 ;       // Twisted Side at phi = 180 deg
  G4VSurface *fSide270 ;       // Twisted Side at phi = 270 deg

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

//---------------------
// inline functions
//---------------------

inline
void G4TwistedTrapezoid::SetFields(G4double phitwist, G4double fHalfZ ,
				   G4double fHalfSideX, G4double fHalfSideY)
{
   fPhiTwist     = phitwist;

   fHalfSides[0]  = fHalfSideX ;
   fHalfSides[1]  = fHalfSideY ;

   fZHalfLength = fHalfZ ;

   //   G4double parity         = (fPhiTwist > 0 ? 1 : -1); 

}

#endif
