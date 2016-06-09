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
// $Id: G4FlatSurface.hh,v 1.4 2004/05/24 12:09:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4FlatSurface
//
// Class description:
//
//  Class describing a flat boundary surface for G4VSolid.

// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------
#ifndef __G4FLATSURFACE__
#define __G4FLATSURFACE__

#include "G4VSurface.hh"

// class G4TwistedTubs;

class G4FlatSurface : public G4VSurface
{
  public:  // with description

   G4FlatSurface(const G4String         &name,
                 const G4RotationMatrix &rot,
                 const G4ThreeVector    &tlate,
                 const G4ThreeVector    &n,
                 const EAxis             axis1 = kRho, // RHO axis !
                 const EAxis             axis2 = kPhi, // PHI axis !
                       G4double          axis0min = -kInfinity,
                       G4double          axis1min = -kInfinity,
                       G4double          axis0max = kInfinity,
                       G4double          axis1max = kInfinity );
                       
   G4FlatSurface( const G4String        &name,
                        G4double         EndInnerRadius[2],
                        G4double         EndOuterRadius[2],
                        G4double         DPhi,
                        G4double         EndPhi[2],
                        G4double         EndZ[2], 
                        G4int            handedness ) ;

   virtual ~G4FlatSurface();
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
                                                  
  protected:  // with description

   virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                   G4bool withTol = true) ;

  private:

   virtual void SetCorners();
   virtual void SetBoundaries();
   
};

#endif
