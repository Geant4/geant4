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
// $Id: G4TwistedBoxSide.hh,v 1.4 2004/12/08 10:20:33 link Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrapSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a box.

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDBOXSIDE__
#define __G4TWISTEDBOXSIDE__

#include "G4VSurface.hh"

#include <vector>

class G4TwistedBoxSide : public G4VSurface
{
  public:  // with description
   
    G4TwistedBoxSide(const G4String     &name,
                           G4double      PhiTwist,
                           G4double      fDx,
                           G4double      fDy,
                           G4double      fDz,
                           G4double      AngleSide);
  
    virtual ~G4TwistedBoxSide();
   
    virtual G4ThreeVector  GetNormal(const G4ThreeVector &xx,
                                            G4bool isGlobal = false) ;   
   
    virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                    const G4ThreeVector &gv,
                                          G4ThreeVector  gxx[],
                                          G4double  distance[],
                                          G4int     areacode[],
                                          G4bool    isvalid[],
                                          EValidate validate=kValidateWithTol);
                                                  
    virtual G4int DistanceToSurface(const G4ThreeVector &gp,
                                          G4ThreeVector  gxx[],
                                          G4double       distance[],
                                          G4int          areacode[]);

  private:

    virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                    G4bool         withTol = true);
    virtual void SetCorners();
    virtual void SetBoundaries();

    void GetPhiUAtX( G4ThreeVector p, G4double &phi, G4double &u);
    G4ThreeVector ProjectPoint(const G4ThreeVector &p,
                                     G4bool isglobal = false);

    inline G4ThreeVector SurfacePoint( G4double phi, G4double u);
    inline G4ThreeVector NormAng( G4double phi, G4double u );

  private:

    G4double       fDx ;
    G4double       fDy ;
    G4double       fDz ;
    G4double       fPhiTwist ;
    G4double       fAngleSide ;
};   

//========================================================
// inline functions
//========================================================

inline
G4ThreeVector G4TwistedBoxSide::SurfacePoint( G4double phi, G4double u ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( fDx * std::cos(phi) - u * std::sin(phi),
                            fDx * std::sin(phi) + u * std::cos(phi),
                            2*fDz*phi/fPhiTwist );
  return SurfPoint ;
}

inline
G4ThreeVector G4TwistedBoxSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface

  G4double L = 2*fDz ;
  G4ThreeVector nvec( L*std::cos(phi), L*std::sin(phi), fPhiTwist*u);

  return nvec.unit() ;
}


#endif
