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
// $Id: G4TwistedTrapSide.hh,v 1.7 2004/12/08 10:20:34 link Exp $
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
//  Class describing a twisted boundary surface for a trapezoid.

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAPSIDE__
#define __G4TWISTEDTRAPSIDE__

#include "G4VSurface.hh"

#include <vector>

class G4TwistedTrapSide : public G4VSurface
{
  public:  // with description
   
    G4TwistedTrapSide(const G4String     &name,
                            G4double      PhiTwist,
                            G4double      pDx1,
                            G4double      pDx2,
                            G4double      pDy,
                            G4double      pDz,
                            G4double      AngleSide);
  
    virtual ~G4TwistedTrapSide();
   
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



  private:

    virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                    G4bool         withTol = true);
    virtual void SetCorners();
    virtual void SetBoundaries();

    void GetPhiUAtX(G4ThreeVector p, G4double &phi, G4double &u);
    G4ThreeVector ProjectPoint(const G4ThreeVector &p,
                                     G4bool isglobal = false);

    inline G4ThreeVector SurfacePoint(G4double phi, G4double u);
    inline G4ThreeVector NormAng(G4double phi, G4double u);
    inline G4double Xcoef(G4double u);    // to calculate the w(u) function

  private:

    G4double fDx1;        // Half-length along x of the side at y=-fDy1
                          // (d in surface equation)
    G4double fDx2;        // Half-length along x of the side at y=+fDy1
                          // (a in surface equation)
    G4double fDy;         // Half-length along y of the face
                          // (b in surface equation)
    G4double fDz;         // Half-length along the z axis
                          // (L in surface equation)
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

};   

//========================================================
// inline functions
//========================================================

inline
G4double G4TwistedTrapSide::Xcoef(G4double u)

{
  // attention a = 2 fDx1 
  //           d = 2 fDx2
  //           b = 2 fDy
  //           L = 2 fDz

  return  fDx2 + ( fDx1-fDx2)/2 - u * (fDx1-fDx2)/(2*fDy) ;
}

inline
G4ThreeVector G4TwistedTrapSide::SurfacePoint( G4double phi, G4double u ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( Xcoef(u) * std::cos(phi) - u * std::sin(phi),
                            Xcoef(u) * std::sin(phi) + u * std::cos(phi),
                            2*fDz*phi/fPhiTwist );
  return SurfPoint ;

}

inline
G4ThreeVector G4TwistedTrapSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface

  G4double L = 2*fDz ;
  G4double a = 2*fDx2 ;
  G4double d = 2*fDx1 ;
  G4double b = 2*fDy ;
  G4double dphi = fPhiTwist ;
  G4double b2 = b*b ;
  G4double d2 = d*d ;
  G4double a2 = a*a ;
  G4ThreeVector nvec( ( L * (2*b*std::cos(phi) + (a-d)*std::sin(phi)))/2.,
                      -(L*((a-d)*std::cos(phi) - 2*b*std::sin(phi)))/2.,
                      (dphi*(-b*d2 + 8*b2*u - 4*a*d*u + 2*d2*u + a2*(b + 2*u)))
                      /(8.*b) ) ;
  return nvec.unit();
}


#endif
