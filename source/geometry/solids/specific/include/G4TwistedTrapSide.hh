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
// $Id: G4TwistedTrapSide.hh,v 1.1 2004-07-29 15:11:16 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//  Class describing a twisted boundary surface for G4VSolid.

// Author: 
//   Oliver Link
//
// --------------------------------------------------------------------
#ifndef __G4TWISTEDTRAPSIDE__
#define __G4TWISTEDTRAPSIDE__

#include "G4VSurface.hh"

class G4TwistedTrapSide : public G4VSurface
{
  public:  // with description
   
  G4TwistedTrapSide(const G4String     &name,
		    G4double      EndInnerRadius[2],
		    G4double      EndOuterRadius[2],
		    G4double      DPhi,
		    G4double      EndPhi[2],
		    G4double      EndZ[2], 
		    G4double      InnerRadius,
		    G4double      OuterRadius,
		    G4double      Kappa,
		    G4double      PhiTwist,
		    G4double      halfzlen,
		    G4double      HalfSides[2],
		    G4int         handedness);

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
 
   inline G4ThreeVector ProjectAtPXPZ(const G4ThreeVector &p,
                                            G4bool isglobal = false) const ;

 private:

  inline G4ThreeVector NormAng( G4double phi, G4double psi ) ;
  inline G4ThreeVector SurfacePoint( G4double phi, G4double psi) ;
  inline void GetPhiPsiAtX(G4ThreeVector x, G4double &phi, G4double &psi) ;
    
  virtual G4double DistanceToPlane(const G4ThreeVector &p,
				   const G4ThreeVector &A,
				   const G4ThreeVector &B,
				   const G4ThreeVector &C,
				   const G4ThreeVector &D,
				   const G4int          parity,
				   G4ThreeVector &xx,
				   G4ThreeVector &n);
  
   virtual G4int GetAreaCode(const G4ThreeVector &xx, 
                                   G4bool         withTol = true);

   virtual void SetCorners();

   virtual void SetCorners(  G4double      endInnerRad[2],
                             G4double      endOuterRad[2],
                             G4double      endPhi[2],
                             G4double      endZ[2] ) ;

   virtual void SetBoundaries();

  private:

  G4double       fKappa;          // tan(TwistedAngle/2)/HalfLenZ;
  G4double       fSideX;
  G4double       fSideY;
  G4double       fZHalfLength ;
  G4double       fPhiTwist ;

};   


//========================================================
// inline functions
//========================================================

inline 
void G4TwistedTrapSide::GetPhiPsiAtX(G4ThreeVector x, G4double &phi, G4double &psi)
{
  phi = x.z()/(2*fZHalfLength)*fPhiTwist ;
  psi = atan( ( fSideY/2 * cos(phi) - x.x() )/ ( fSideY/2 * sin(phi) ) ) ;
}

inline
G4ThreeVector G4TwistedTrapSide::SurfacePoint( G4double phi, G4double psi) 
{
  G4ThreeVector vec( fSideY/2 * cos(phi) - fSideY/2 * sin(phi)*tan(psi),
		     fSideY/2 * sin(phi) + fSideY/2 * cos(phi)*tan(psi),
		     2*fZHalfLength*phi/fPhiTwist );
  return vec ;
}



inline 
G4ThreeVector G4TwistedTrapSide::NormAng( G4double phi, G4double psi ) 
{
  // function to calculate the norm at a given point on the surface
  G4double L = 2*fZHalfLength ;
  G4ThreeVector nvec( L*cos(phi), L*sin(phi), fSideY/2*fPhiTwist*tan(psi));
  return nvec / nvec.mag() ;
}


inline
G4ThreeVector G4TwistedTrapSide::ProjectAtPXPZ(const G4ThreeVector &p, 
                                                    G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }
  G4ThreeVector xx(p.x(), p.x() * fKappa * p.z(), p.z());
  if (isglobal) {
     return (fRot * xx + fTrans);
  } else {
     return xx;
  }
}


#endif
