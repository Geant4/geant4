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
// $Id: G4FlatTrapSide.hh,v 1.3 2004/11/13 18:26:24 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4FlatTrapSurface
//
// Class description:
//
//  Class describing a flat boundary surface for a trapezoid.

// Author:
//
//   27-Oct-2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4FLATTRAPSIDE__
#define __G4FLATTRAPSIDE__

#include "G4VSurface.hh"

class G4FlatTrapSide : public G4VSurface
{
  public:  // with description

   G4FlatTrapSide( const G4String& name,
                         G4double  PhiTwist,
                         G4double  pDx1,
                         G4double  pDx2,
                         G4double  pDy,
                         G4double  pDz,
                         G4int     handedness  );
   virtual ~G4FlatTrapSide();

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

  private:
  
   G4double fDx1 ;
   G4double fDx2 ;
   G4double fDy ;
   G4double fDz ;
   G4double fPhiTwist ;
};

#endif
