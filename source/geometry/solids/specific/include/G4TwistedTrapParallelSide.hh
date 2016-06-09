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
// $Id: G4TwistedTrapParallelSide.hh,v 1.3 2005/03/18 17:11:53 gcosmo Exp $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistedTrapParallelSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a trapezoid.

// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAPPARALLELSIDE__
#define __G4TWISTEDTRAPPARALLELSIDE__

#include "G4VSurface.hh"

#include <vector>

class G4TwistedTrapParallelSide : public G4VSurface
{
  public:  // with description
   
    G4TwistedTrapParallelSide(const G4String     &name,
			   G4double      PhiTwist,    // twist angle
			   G4double      pDz,         // half z lenght
			   G4double      pTheta,      // direction between end planes
			   G4double      pPhi,        // defined by polar and azimutal angles.
			   G4double      pDy1,        // half y length at -pDz
			   G4double      pDx1,        // half x length at -pDz,-pDy
			   G4double      pDx2,        // half x length at -pDz,+pDy
			   G4double      pDy2,        // half y length at +pDz
			   G4double      pDx3,        // half x length at +pDz,-pDy
			   G4double      pDx4,        // half x length at +pDz,+pDy
			   G4double      pAlph,       // tilt angle at +pDz
                           G4double      AngleSide    // parity
			   );
  
    virtual ~G4TwistedTrapParallelSide();
   
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
    inline G4double GetValueB(G4double phi) ;
    inline G4double GetUmin(G4double phi) ;
    inline G4double GetUmax(G4double phi) ;

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
    G4double fTAlph ;    // std::tan(fAlph)
    
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

    G4double fa1md1 ;   // 2 fDx2 - 2 fDx1  == a1 - d1
    G4double fa2md2 ;  // 2 fDx4 - 2 fDx3 

    G4double fdeltaX ;
    G4double fdeltaY ;

    G4double fDy2plus1 ;  // fDy2 + fDy1  == b2/2 + b1/2
    G4double fDy2minus1 ; // fDy2 - fDy1          -

  // temporary
  G4double fDy ;

};   

//========================================================
// inline functions
//========================================================

inline 
G4double G4TwistedTrapParallelSide::GetUmin(G4double phi)
{
  return ( - (2*fDx2*(fPhiTwist - 2*phi) + 2*fDx4*(fPhiTwist + 2*phi) + (2*fDy1*(fPhiTwist - 2*phi) + 2*fDy2*(fPhiTwist + 2*phi))*fTAlph)/(4.*fPhiTwist) ) ;
}

inline
G4double G4TwistedTrapParallelSide::GetUmax(G4double phi)
{
  return   (2*fDx2*(fPhiTwist - 2*phi) + 2*fDx4*(fPhiTwist + 2*phi) - (2*fDy1*(fPhiTwist - 2*phi) + 2*fDy2*(fPhiTwist + 2*phi))*fTAlph)/(4.*fPhiTwist) ;
}

inline 
G4double G4TwistedTrapParallelSide::GetValueB(G4double phi) 
{
  return ( fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ) ;
}


inline
G4double G4TwistedTrapParallelSide::Xcoef(G4double phi)
{
  
  return GetValueB(phi)/2. ;

}

inline
G4ThreeVector G4TwistedTrapParallelSide::SurfacePoint( G4double phi, G4double u ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( ( u + fdeltaX*phi/fPhiTwist ) * std::cos(phi) 
			    - ( Xcoef(phi) + fdeltaY*phi/fPhiTwist ) * std::sin(phi),
                            ( u + fdeltaX*phi/fPhiTwist ) * std::sin(phi)
			    + ( Xcoef(phi) + fdeltaY*phi/fPhiTwist )  * std::cos(phi),
			    2*fDz*phi/fPhiTwist  );
  return SurfPoint ;

}

inline
G4ThreeVector G4TwistedTrapParallelSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface
  // replace a1-d1

  G4ThreeVector nvec( 
		     -2*fDz*std::sin(phi),
		     2*fDz*std::cos(phi),
		     fDy1 - fDy2 - fdeltaY - fdeltaX*phi - fPhiTwist*u
 ) ;


  return nvec.unit();
}


#endif
