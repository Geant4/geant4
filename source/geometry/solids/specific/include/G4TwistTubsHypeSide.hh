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
// $Id: G4TwistTubsHypeSide.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistTubsHypeSide
//
// Class description:
//
//  Class describing a hyperbolic boundary surface for a cylinder.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4TWISTTUBSHYPESIDE__
#define __G4TWISTTUBSHYPESIDE__

#include "G4VTwistSurface.hh"
#include "G4Integrator.hh"
#include "G4SimpleIntegration.hh"

class G4TwistTubsHypeSide : public G4VTwistSurface
{
  public:  // with description
                       
   G4TwistTubsHypeSide(const G4String         &name,
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
                             
  G4TwistTubsHypeSide(const G4String  &name,
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

   virtual ~G4TwistTubsHypeSide();

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
   
   virtual G4double GetRhoAtPZ(const G4ThreeVector &p,
                                     G4bool isglobal = false) const ;
   
   virtual G4ThreeVector SurfacePoint(G4double, G4double,
                                      G4bool isGlobal = false) ;  
   virtual G4double GetBoundaryMin(G4double phi) ;
   virtual G4double GetBoundaryMax(G4double phi) ;
   virtual G4double GetSurfaceArea() ;
   virtual void GetFacets( G4int m, G4int n, G4double xyz[][3],
                           G4int faces[][4], G4int iside ) ;
 
  public:  // without description

   G4TwistTubsHypeSide(__void__&);
     // Fake default constructor for usage restricted to direct object
     // persistency for clients requiring preallocation of memory for
     // persistifiable objects.

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
   
   G4double          fKappa;        // std::tan(TwistedAngle/2)/HalfLenZ;
   G4double          fTanStereo;    // std::tan(StereoAngle)
   G4double          fTan2Stereo;   // std::tan(StereoAngle)**2
   G4double          fR0;           // radius at z = 0
   G4double          fR02;          // radius**2 at z = 0
   G4double          fDPhi ;        // segment

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
G4double G4TwistTubsHypeSide::GetRhoAtPZ(const G4ThreeVector &p,
                                               G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }
  return std::sqrt(fR02 + tmpp.z() * tmpp.z() * fTan2Stereo); 
}

inline
G4ThreeVector G4TwistTubsHypeSide::
SurfacePoint(G4double phi , G4double z , G4bool isGlobal)
{
  G4double rho = std::sqrt(fR02 + z * z * fTan2Stereo) ;

  G4ThreeVector SurfPoint (rho*std::cos(phi), rho*std::sin(phi), z) ;

  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTubsHypeSide::GetBoundaryMin(G4double z)
{
  G4ThreeVector ptmp(0,0,z) ;  // temporary point with z Komponent only
  G4ThreeVector lowerlimit;    // lower phi-boundary limit at z = ptmp.z()
  lowerlimit = GetBoundaryAtPZ(sAxis0 & sAxisMin, ptmp);
  return  std::atan2( lowerlimit.y(), lowerlimit.x() ) ;  
}

inline
G4double G4TwistTubsHypeSide::GetBoundaryMax(G4double z )
{
  G4ThreeVector ptmp(0,0,z) ;  // temporary point with z Komponent only
  G4ThreeVector upperlimit;    // upper phi-boundary limit at z = ptmp.z()
  upperlimit = GetBoundaryAtPZ(sAxis0 & sAxisMax, ptmp);
  return   std::atan2( upperlimit.y(), upperlimit.x() ) ;
}

inline
G4double G4TwistTubsHypeSide::GetSurfaceArea()
{
  // approximation with tube surface

  return ( fAxisMax[1] - fAxisMin[1] ) * fR0 * fDPhi ;
}

#endif
