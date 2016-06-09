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
// $Id: G4HyperbolicSurface.hh,v 1.4 2004/05/24 12:09:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4HyperbolicSurface
//
// Class description:
//
//  Class describing a hyperbolic boundary surface for G4VSolid.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4HYPERBOLICSURFACE__
#define __G4HYPERBOLICSURFACE__

#include "G4VSurface.hh"

class G4HyperbolicSurface : public G4VSurface
{
  public:  // with description
                       
   G4HyperbolicSurface(const G4String         &name,
                       const G4RotationMatrix &rot,  // 0.5*(phi-width segment)
                       const G4ThreeVector    &tlate,
                       const G4int     handedness,// R-hand = 1, L-hand = -1
                       const G4double  kappa,     // tan(TwistAngle/2)/fZHalfLen
                       const G4double  tanstereo, // tan(stereo angle)
                       const G4double  r0,        // radius at z = 0
                       const EAxis     axis0 = kPhi,
                       const EAxis     axis1 = kZAxis,
                             G4double  axis0min = -kInfinity,
                             G4double  axis1min = -kInfinity,
                             G4double  axis0max = kInfinity,
                             G4double  axis1max = kInfinity); 
                             
  G4HyperbolicSurface(const G4String  &name,
                            G4double   EndInnerRadius[2],
                            G4double   EndOuterRadius[2],
                            G4double   DPhi,
                            G4double   EndPhi[2],
                            G4double   EndZ[2], 
                            G4double   InnerRadius,
                            G4double   OuterRadius,
                            G4double   Kappa,
                            G4double   TanInnerStereo,
                            G4double   TanOuterStereo,
                            G4int      handedness) ;

   virtual ~G4HyperbolicSurface();

   virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                   const G4ThreeVector &gv,
                                         G4ThreeVector  gxx[],
                                         G4double       distance[],
                                         G4int          areacode[],
                                         G4bool         isvalid[],
                                   EValidate validate = kValidateWithTol);
                                                   
   virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                         G4ThreeVector  gxx[],
                                         G4double       distance[],
                                         G4int          areacode[]);
 
   virtual G4ThreeVector GetNormal(const G4ThreeVector &xx,
                                         G4bool isGlobal = false) ;
   virtual EInside Inside(const G4ThreeVector &gp) ;
   
   inline virtual G4double GetRhoAtPZ(const G4ThreeVector &p,
                                            G4bool isglobal = false) const ;
   
  private:

   virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                   G4bool withTol = true);
   virtual G4int GetAreaCodeInPhi(const G4ThreeVector &xx, 
                                        G4bool withTol = true);
   virtual void SetCorners();

   virtual void SetCorners(G4double         EndInnerRadius[2],
                           G4double         EndOuterRadius[2],
                           G4double         DPhi,
                           G4double         EndPhi[2],
                           G4double         EndZ[2]);
   virtual void SetBoundaries();

  private:
   
   G4double          fKappa;        // tan(TwistedAngle/2)/HalfLenZ;
   G4double          fTanStereo;    // tan(StereoAngle)
   G4double          fTan2Stereo;   // tan(StereoAngle)**2
   G4double          fR0;           // radius at z = 0
   G4double          fR02;          // radius**2 at z = 0
   class Insidetype
   {
     public:
       G4ThreeVector gp;
       EInside       inside;
   };
   Insidetype fInside;
};

//========================================================
// inline functions
//========================================================

inline
G4double G4HyperbolicSurface::GetRhoAtPZ(const G4ThreeVector &p,
                                               G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }
  return sqrt(fR02 + tmpp.z() * tmpp.z() * fTan2Stereo); 
}

#endif
