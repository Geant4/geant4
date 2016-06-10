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
// $Id: G4TwistTubsFlatSide.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistTubsFlatSide
//
// Class description:
//
//  Class describing a flat boundary surface for a cylinder.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4TWISTTUBSFLATSIDE__
#define __G4TWISTTUBSFLATSIDE__

#include "G4VTwistSurface.hh"

class G4TwistTubsFlatSide : public G4VTwistSurface
{
  public:  // with description

   G4TwistTubsFlatSide(const G4String         &name,
                 const G4RotationMatrix &rot,
                 const G4ThreeVector    &tlate,
                 const G4ThreeVector    &n,
                 const EAxis             axis1 = kRho, // RHO axis !
                 const EAxis             axis2 = kPhi, // PHI axis !
                       G4double          axis0min = -kInfinity,
                       G4double          axis1min = -kInfinity,
                       G4double          axis0max = kInfinity,
                       G4double          axis1max = kInfinity );
                       
   G4TwistTubsFlatSide( const G4String        &name,
                        G4double         EndInnerRadius[2],
                        G4double         EndOuterRadius[2],
                        G4double         DPhi,
                        G4double         EndPhi[2],
                        G4double         EndZ[2], 
                        G4int            handedness ) ;

   virtual ~G4TwistTubsFlatSide();
   virtual G4ThreeVector  GetNormal(const G4ThreeVector & /* xx */ ,
                                          G4bool isGlobal = false);
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
                                                  
  virtual G4ThreeVector SurfacePoint(G4double, G4double,
                                     G4bool isGlobal = false ) ;  
  virtual G4double GetBoundaryMin(G4double phi) ;
  virtual G4double GetBoundaryMax(G4double phi) ;
  virtual G4double GetSurfaceArea() { return fSurfaceArea ; }
  virtual void GetFacets( G4int m, G4int n, G4double xyz[][3],
                          G4int faces[][4], G4int iside ) ;

  G4double fSurfaceArea ;

  public:  // without description

   G4TwistTubsFlatSide(__void__&);
     // Fake default constructor for usage restricted to direct object
     // persistency for clients requiring preallocation of memory for
     // persistifiable objects.

  protected:  // with description

   virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                   G4bool withTol = true) ;

  private:

   virtual void SetCorners();
   virtual void SetBoundaries();
   
};

inline G4ThreeVector G4TwistTubsFlatSide::
SurfacePoint(G4double phi , G4double rho , G4bool isGlobal )
{
  G4ThreeVector SurfPoint (rho*std::cos(phi) , rho*std::sin(phi) , 0);

  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}

inline
G4double G4TwistTubsFlatSide::GetBoundaryMin(G4double)
{
  G4ThreeVector dphimin = GetCorner(sC0Max1Min);
  return  std::atan2( dphimin.y(), dphimin.x() );  
}

inline
G4double G4TwistTubsFlatSide::GetBoundaryMax(G4double)
{
  G4ThreeVector dphimax = GetCorner(sC0Max1Max);   
  return  std::atan2( dphimax.y(), dphimax.x() );  
}

#endif
