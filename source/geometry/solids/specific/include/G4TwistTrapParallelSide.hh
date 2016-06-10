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
// $Id: G4TwistTrapParallelSide.hh 67973 2013-03-13 10:16:25Z gcosmo $
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4TwistTrapParallelSide
//
// Class description:
//
//  Class describing a twisted boundary surface for a trapezoid.

// Author:
//
//   27-Oct-2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------
#ifndef __G4TWISTTRAPPARALLELSIDE__
#define __G4TWISTTRAPPARALLELSIDE__

#include "G4VTwistSurface.hh"

#include <vector>

class G4TwistTrapParallelSide : public G4VTwistSurface
{
  public:  // with description
   
    G4TwistTrapParallelSide(const G4String &name,
                              G4double  PhiTwist, // twist angle
                              G4double  pDz,      // half z lenght
                              G4double  pTheta, // direction between end planes
                              G4double  pPhi,   // by polar and azimutal angles
                              G4double  pDy1,     // half y length at -pDz
                              G4double  pDx1,     // half x length at -pDz,-pDy
                              G4double  pDx2,     // half x length at -pDz,+pDy
                              G4double  pDy2,     // half y length at +pDz
                              G4double  pDx3,     // half x length at +pDz,-pDy
                              G4double  pDx4,     // half x length at +pDz,+pDy
                              G4double  pAlph,    // tilt angle at +pDz
                              G4double  AngleSide // parity
                            );
  
    virtual ~G4TwistTrapParallelSide();
   
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

  public:  // without description

    G4TwistTrapParallelSide(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  private:

    virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                    G4bool         withTol = true);
    virtual void SetCorners();
    virtual void SetBoundaries();

    void GetPhiUAtX(G4ThreeVector p, G4double &phi, G4double &u);
    G4ThreeVector ProjectPoint(const G4ThreeVector &p,
                                     G4bool isglobal = false);

    virtual G4ThreeVector SurfacePoint(G4double phi, G4double u,
                                       G4bool isGlobal = false);
    virtual G4double GetBoundaryMin(G4double phi);
    virtual G4double GetBoundaryMax(G4double phi);
    virtual G4double GetSurfaceArea();
    virtual void GetFacets( G4int m, G4int n, G4double xyz[][3],
                            G4int faces[][4], G4int iside );

    inline G4ThreeVector NormAng(G4double phi, G4double u);
    inline G4double GetValueB(G4double phi) ;
    inline G4double Xcoef(G4double u);
      // To calculate the w(u) function

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

    G4double fAlph;
    G4double fTAlph;      // std::tan(fAlph)
    
    G4double fPhiTwist;   // twist angle ( dphi in surface equation)

    G4double fAngleSide;

    G4double fdeltaX;
    G4double fdeltaY;

    G4double fDx4plus2;  // fDx4 + fDx2  == a2/2 + a1/2
    G4double fDx4minus2; // fDx4 - fDx2          -
    G4double fDx3plus1;  // fDx3 + fDx1  == d2/2 + d1/2
    G4double fDx3minus1; // fDx3 - fDx1          -
    G4double fDy2plus1;  // fDy2 + fDy1  == b2/2 + b1/2
    G4double fDy2minus1; // fDy2 - fDy1          -
    G4double fa1md1;     // 2 fDx2 - 2 fDx1  == a1 - d1
    G4double fa2md2;     // 2 fDx4 - 2 fDx3 
};   

//========================================================
// inline functions
//========================================================


inline 
G4double G4TwistTrapParallelSide::GetValueB(G4double phi) 
{
  return ( fDy2plus1 + fDy2minus1 * ( 2 * phi ) / fPhiTwist ) ;
}

inline
G4double G4TwistTrapParallelSide::Xcoef(G4double phi)
{
  return GetValueB(phi)/2. ;
}

inline G4ThreeVector
G4TwistTrapParallelSide::
SurfacePoint( G4double phi, G4double u, G4bool isGlobal ) 
{
  // function to calculate a point on the surface, given by parameters phi,u

  G4ThreeVector SurfPoint ( u*std::cos(phi) - Xcoef(phi)*std::sin(phi)
                          + fdeltaX*phi/fPhiTwist,
                            u*std::sin(phi) + Xcoef(phi)*std::cos(phi)
                          + fdeltaY*phi/fPhiTwist,
                            2*fDz*phi/fPhiTwist  );
  if (isGlobal) { return (fRot * SurfPoint + fTrans); }
  return SurfPoint;
}


inline
G4double G4TwistTrapParallelSide::GetBoundaryMin(G4double phi)
{
  return  -(fPhiTwist*(fDx2 + fDx4 - fDy2plus1*fTAlph)
          + 2*fDx4minus2*phi - 2*fDy2minus1*fTAlph*phi)/(2.*fPhiTwist) ;
}

inline
G4double G4TwistTrapParallelSide::GetBoundaryMax(G4double phi)
{
  return (fDx2 + fDx4 + fDy2plus1*fTAlph)/ 2.
       + ((fDx4minus2 + fDy2minus1*fTAlph)*phi)/fPhiTwist ;
}

inline
G4double  G4TwistTrapParallelSide::GetSurfaceArea()
{
  return 2*fDx4plus2*fDz ;
}

inline
G4ThreeVector G4TwistTrapParallelSide::NormAng( G4double phi, G4double u ) 
{
  // function to calculate the norm at a given point on the surface
  // replace a1-d1

  G4ThreeVector nvec(-2*fDz*std::sin(phi) ,
                      2*fDz*std::cos(phi) ,
                     -(fDy2minus1 + fPhiTwist*u + fdeltaY*std::cos(phi)
                     -fdeltaX*std::sin(phi))) ;
  return nvec.unit();
}

#endif
