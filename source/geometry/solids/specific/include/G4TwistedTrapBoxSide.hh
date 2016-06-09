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
// $Id: G4TwistedTrapBoxSide.hh,v 1.2 2005/03/18 17:11:53 gcosmo Exp $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrapBoxSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a trapezoid.

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAPBOXSIDE__
#define __G4TWISTEDTRAPBOXASIDE__

#include "G4VSurface.hh"

#include <vector>

class G4TwistedTrapBoxSide : public G4VSurface
{
  public:  // with description
   
    G4TwistedTrapBoxSide(const G4String     &name,
			   G4double      PhiTwist,    // twist angle
			   G4double      pDz,         // half z lenght
			   G4double      pTheta,      // direction between end planes
			   G4double      pPhi,        // defined by polar and azimutal angles.
			   G4double      pDy1,        // half y length at -pDz
			   G4double      pDx1,        // half x length at -pDz,-pDy
			   G4double      pDy2,        // half y length at +pDz
			   G4double      pDx2,        // half x length at +pDz,-pDy
			   G4double      pAlph,       // tilt angle at +pDz
                           G4double      AngleSide    // parity
			   );
  
    virtual ~G4TwistedTrapBoxSide();
   
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
    inline G4double Xcoef(G4double u,G4double phi);    // to calculate the w(u) function
    inline G4double GetValueA(G4double phi) ;
    inline G4double GetValueB(G4double phi) ;

  private:

    G4double fTheta;   
    G4double fPhi ;

    G4double fDy1;   
    G4double fDy2;   

    G4double fDx1;     
    G4double fDx2;     

    G4double fDz;         // Half-length along the z axis

    G4double fAlph ;
    G4double fTAlph ;    // std::tan(fAlph)
    
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

    G4double fdeltaX ;
    G4double fdeltaY ;

};   

//========================================================
// inline functions
//========================================================

inline
G4double G4TwistedTrapBoxSide::GetValueA(G4double phi)
{
  return ( fDx1 + fDx2 - (2*(fDx1 - fDx2)*phi)/fPhiTwist ) ;
}

inline 
G4double G4TwistedTrapBoxSide::GetValueB(G4double phi) 
{
  return ( fDy1 + fDy2 - (2*(fDy1 - fDy2)*phi)/fPhiTwist ) ;
}


inline
G4double G4TwistedTrapBoxSide::Xcoef(G4double u, G4double phi)
{
  
  return GetValueA(phi)/2. - u*fTAlph    ;

}

inline
G4ThreeVector G4TwistedTrapBoxSide::SurfacePoint( G4double phi, G4double u ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( ( Xcoef(u,phi) + fdeltaX*phi/fPhiTwist ) * std::cos(phi) 
			    - ( u + fdeltaY*phi/fPhiTwist )          * std::sin(phi),
                            ( Xcoef(u,phi) + fdeltaX*phi/fPhiTwist ) * std::sin(phi)
			    + ( u + fdeltaY*phi/fPhiTwist )          * std::cos(phi),
			    2*fDz*phi/fPhiTwist  );
  return SurfPoint ;

}

inline
G4ThreeVector G4TwistedTrapBoxSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface
  // replace a1-d1

  G4ThreeVector 
    nvec(
	 8*fDz*(std::cos(phi) - fTAlph*std::sin(phi)),
	 8*fDz*(fTAlph*std::cos(phi) + std::sin(phi)),
	 2*fDx1*(2 - fTAlph*(fPhiTwist - 2*phi)) - 2*fDx2*(2 + fTAlph*(fPhiTwist + 2*phi)) 
	 - 4*(fdeltaX + fdeltaY*(fTAlph - phi) + fdeltaX*fTAlph*phi) + 4*fPhiTwist*(1 + fTAlph*fTAlph)*u) ;

  return nvec.unit();
}


#endif
