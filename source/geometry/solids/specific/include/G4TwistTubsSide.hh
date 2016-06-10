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
// $Id: G4TwistTubsSide.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistTubsSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a cylinder.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4TWISTTUBSSIDE__
#define __G4TWISTTUBSSIDE__

#include "G4VTwistSurface.hh"

class G4TwistTubsSide : public G4VTwistSurface
{
  public:  // with description
   
   G4TwistTubsSide(const G4String         &name,
                   const G4RotationMatrix &rot,   // 0.5*(phi-width segment)
                   const G4ThreeVector    &tlate,
                         G4int    handedness, // R-hand = 1, L-hand = -1
                   const G4double kappa,      // tan(TwistAngle/2)/fZHalfLen
                   const EAxis    axis0 = kXAxis,
                   const EAxis    axis1 = kZAxis,
                         G4double axis0min = -kInfinity,
                         G4double axis1min = -kInfinity,
                         G4double axis0max = kInfinity,
                         G4double axis1max = kInfinity );
    
   G4TwistTubsSide(const G4String     &name,
                         G4double      EndInnerRadius[2],
                         G4double      EndOuterRadius[2],
                         G4double      DPhi,
                         G4double      EndPhi[2],
                         G4double      EndZ[2], 
                         G4double      InnerRadius,
                         G4double      OuterRadius,
                         G4double      Kappa,
                         G4int         handedness);

   virtual ~G4TwistTubsSide();
   
   virtual G4ThreeVector  GetNormal(const G4ThreeVector &xx,
                                          G4bool isGlobal = false) ;   
   
   virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                   const G4ThreeVector &gv,
                                         G4ThreeVector  gxx[],
                                         G4double  distance[],
                                         G4int     areacode[],
                                         G4bool    isvalid[],
                                         EValidate validate = kValidateWithTol);
                                                  
   virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                         G4ThreeVector  gxx[],
                                         G4double       distance[],
                                         G4int          areacode[]);
 
   inline G4ThreeVector ProjectAtPXPZ(const G4ThreeVector &p,
                                            G4bool isglobal = false) const ;

   virtual G4ThreeVector SurfacePoint(G4double, G4double,
                                      G4bool isGlobal = false) ;  
   virtual G4double GetBoundaryMin(G4double phi) ;
   virtual G4double GetBoundaryMax(G4double phi) ;
   virtual G4double GetSurfaceArea() ;
   virtual void GetFacets( G4int m, G4int n, G4double xyz[][3],
                           G4int faces[][4], G4int iside ) ;

 public:  // without description

   G4TwistTubsSide(__void__&);
     // Fake default constructor for usage restricted to direct object
     // persistency for clients requiring preallocation of memory for
     // persistifiable objects.

 private:

   virtual G4double DistanceToPlane(const G4ThreeVector &p,
                                    const G4ThreeVector &A,
                                    const G4ThreeVector &B,
                                    const G4ThreeVector &C,
                                    const G4ThreeVector &D,
                                    const G4int          parity,
                                          G4ThreeVector &xx,
                                          G4ThreeVector &n);

   virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                   G4bool         withTol = true);

   virtual void SetCorners();

   virtual void SetCorners(  G4double      endInnerRad[2],
                             G4double      endOuterRad[2],
                             G4double      endPhi[2],
                             G4double      endZ[2] ) ;

   virtual void SetBoundaries();

  private:

   G4double       fKappa;          // std::tan(TwistedAngle/2)/HalfLenZ;
};   


//========================================================
// inline functions
//========================================================

inline
G4ThreeVector G4TwistTubsSide::ProjectAtPXPZ(const G4ThreeVector &p, 
                                                   G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }
  G4ThreeVector xx(p.x(), p.x() * fKappa * p.z(), p.z());
  if (isglobal) { return (fRot * xx + fTrans); }
  return xx;
}

inline
G4ThreeVector
G4TwistTubsSide::SurfacePoint(G4double x, G4double z, G4bool isGlobal)
{
  G4ThreeVector SurfPoint( x , x * fKappa * z , z ) ;

  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTubsSide::GetBoundaryMin(G4double)
{
  return  fAxisMin[0] ;  // inner radius at z = 0
}

inline
G4double G4TwistTubsSide::GetBoundaryMax(G4double)
{
  return  fAxisMax[0] ;  // outer radius at z = 0
}

inline
G4double G4TwistTubsSide::GetSurfaceArea()
{
  // approximation only
  return ( fAxisMax[0] - fAxisMin[0] ) * ( fAxisMax[1] - fAxisMin[1] ) ;
}

#endif
