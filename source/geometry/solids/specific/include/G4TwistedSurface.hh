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
// $Id: G4TwistedSurface.hh,v 1.4 2004/05/24 12:09:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedSurface
//
// Class description:
//
//  Class describing a twisted boundary surface for G4VSolid.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4TWISTEDSURFACE__
#define __G4TWISTEDSURFACE__

#include "G4VSurface.hh"

class G4TwistedSurface : public G4VSurface
{
  public:  // with description
   
   G4TwistedSurface(const G4String         &name,
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
    
   G4TwistedSurface(const G4String     &name,
                          G4double      EndInnerRadius[2],
                          G4double      EndOuterRadius[2],
                          G4double      DPhi,
                          G4double      EndPhi[2],
                          G4double      EndZ[2], 
                          G4double      InnerRadius,
                          G4double      OuterRadius,
                          G4double      Kappa,
                          G4int         handedness);

   virtual ~G4TwistedSurface();
   
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

   G4double       fKappa;          // tan(TwistedAngle/2)/HalfLenZ;
};   


//========================================================
// inline functions
//========================================================

inline
G4ThreeVector G4TwistedSurface::ProjectAtPXPZ(const G4ThreeVector &p, 
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
  if (isglobal) {
     return (fRot * xx + fTrans);
  } else {
     return xx;
  }
}

#endif
