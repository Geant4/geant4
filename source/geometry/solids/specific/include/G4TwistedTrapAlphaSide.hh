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
// $Id: G4TwistedTrapAlphaSide.hh,v 1.1 2005-02-15 13:48:24 link Exp $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrapAlphaSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a trapezoid.

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAPALPHASIDE__
#define __G4TWISTEDTRAPALPHASIDE__

#include "G4VSurface.hh"

#include <vector>

class G4TwistedTrapAlphaSide : public G4VSurface
{
  public:  // with description
   
    G4TwistedTrapAlphaSide(const G4String     &name,
			   G4double      PhiTwist,    // twist angle
			   G4double      pDz,         // half z lenght
			   G4double      pTheta,      // direction between end planes
			   G4double      pPhi,        // defined by polar and azimutal angles.
			   G4double      pDy1,        // half y length at -pDz
			   G4double      pDx1,        // half x length at -pDz,-pDy
			   G4double      pDx2,        // half x length at -pDz,+pDy
			   G4double      pAlph1,      // tilt angle at -pDz
			   G4double      pDy2,        // half y length at +pDz
			   G4double      pDx3,        // half x length at +pDz,-pDy
			   G4double      pDx4,        // half x length at +pDz,+pDy
			   G4double      pAlph2,      // tilt angle at +pDz
                           G4double      AngleSide    // parity
			   );
  
    virtual ~G4TwistedTrapAlphaSide();
   
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
    inline G4double GetValueD(G4double phi) ;

  private:

    G4double fTheta;   
    G4double fPhi ;

    G4double fDy1;   
    G4double fDx1;     
    G4double fDx2;     

    G4double fDy2;   
    G4double fDx3;     
    G4double fDx4;     

    G4double fDz;         // Half-length along the z axis

    G4double fAlph ;
    G4double fTAlph ;    // tan(fAlph)
    
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

    G4double fDx4plus2 ;  // fDx4 + fDx2  == a2/2 + a1/2
    G4double fDx4minus2 ; // fDx4 - fDx2          -
    G4double fDx3plus1  ; // fDx3 + fDx1  == d2/2 + d1/2
    G4double fDx3minus1 ; // fDx3 - fDx1          -
    G4double fDy2plus1 ;  // fDy2 + fDy1  == b2/2 + b1/2
    G4double fDy2minus1 ; // fDy2 - fDy1          -
    G4double fa1md1 ;   // 2 fDx2 - 2 fDx3  == a1 - d1
 


    G4double fdeltaX ;
    G4double fdeltaY ;


  // temporary
  G4double fDy ;

};   

//========================================================
// inline functions
//========================================================

inline
G4double G4TwistedTrapAlphaSide::GetValueA(G4double phi)
{
  return ( fDx4plus2 + fDx4minus2 * ( 2 * phi ) / fPhiTwist  ) ;
}

inline
G4double G4TwistedTrapAlphaSide::GetValueD(G4double phi) 
{
  return ( fDx3plus1 + fDx3minus1 * ( 2 * phi) / fPhiTwist  ) ;
} 

inline 
G4double G4TwistedTrapAlphaSide::GetValueB(G4double phi) 
{
  return ( fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ) ;
}


inline
G4double G4TwistedTrapAlphaSide::Xcoef(G4double u, G4double phi)
{
  
  return GetValueA(phi)/2. + (GetValueD(phi)-GetValueA(phi))/4. 
    - u*( ( GetValueD(phi)-GetValueA(phi) ) / ( 2 * GetValueB(phi) ) + fTAlph )   ;

}

inline
G4ThreeVector G4TwistedTrapAlphaSide::SurfacePoint( G4double phi, G4double u ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( ( Xcoef(u,phi) + fdeltaX*phi ) * std::cos(phi) 
			    - ( u + fdeltaY*phi )          * std::sin(phi),
                            ( Xcoef(u,phi) + fdeltaX*phi ) * std::sin(phi)
			    + ( u + fdeltaY*phi )          * std::cos(phi),
			    2*fDz*phi/fPhiTwist  );
  return SurfPoint ;

}

inline
G4ThreeVector G4TwistedTrapAlphaSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface
  // replace a1-d1

  G4double L = 2*fDz ;

 
  G4ThreeVector nvec( 8*b1*L*(2*b1*(std::cos(phi)-std::sin(phi)*fTAlph) + fa1md1*std::sin(phi)),
		      8*b1*L*(2*b1*std::sin(phi) + std::cos(phi)*(-fa1md1 + 2*b1*fTAlph)),
		      b1*(4*(a1-a2)*b1 + 4*b1*(d1-d2) + a1*(a1+a2)*dphi -
			  (a2+d1)*d1*dphi + fa1md1*d2*dphi  - 
			  16*b1*dx + 8*fa1md1*dy  + 
			  2*(-(fa1md1*(a1 - a2 + d1 - d2 - 4*dx)) + 8*b1*dy)*phi) + 
		      4*(4*b1*b1 + (fa1md1)*(fa1md1))*dphi*u - 
		      2*b1*fTAlph*(b1*((a1 + a2 + d1 + d2)*dphi + 8*dy - 
				       2*(a1 - a2 + d1 - d2 - 4*dx)*phi) + 8*(fa1md1)*dphi*u - 
				   8*b1*dphi*u*fTAlph) ) ;


  return nvec.unit();
}


#endif
