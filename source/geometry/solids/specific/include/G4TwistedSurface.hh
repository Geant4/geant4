// $Id: G4TwistedSurface.hh,v 1.1 2003-11-12 15:43:00 link Exp $
#ifndef __G4TWISTEDSURFACE__
#define __G4TWISTEDSURFACE__
//*************************************************************************
//* --------------------
//* G4TwistedSurface
//* --------------------
//* (Description)
//* 	Class for twisted boundary surface for J4Solid.
//*     
//* (Update Record)
//*	2002/08/01  K.Hoshina	Original version.
//*************************************************************************


#include "G4VSurface.hh"

class G4TwistedTubs;

class G4TwistedSurface : public G4VSurface
{
public:
   
   G4TwistedSurface(const G4String         &name,
                    const G4RotationMatrix &rot,        // rotation of 0.5*(phi-width of a segment)
                    const G4ThreeVector    &tlate,
                    G4int                   handedness, // right hand = 1, left hand = -1
                    const G4double          kappa,      // tan(TwistAngle/2)/fZHalfLen
                    const EAxis             axis0 = kXAxis,
                    const EAxis             axis1 = kZAxis,
                          G4double          axis0min = -kInfinity,
                          G4double          axis1min = -kInfinity,
                          G4double          axis0max = kInfinity,
                          G4double          axis1max = kInfinity );
                          
   G4TwistedSurface(const G4String &name, G4TwistedTubs *solid, G4int handedness);

   virtual ~G4TwistedSurface();
   
   virtual G4ThreeVector  GetNormal(const G4ThreeVector &xx, G4bool isGlobal = FALSE) ;   
   
   virtual G4int          DistanceToSurface(const G4ThreeVector &gp,
                                            const G4ThreeVector &gv,
                                                  G4ThreeVector  gxx[],
                                                  G4double       distance[],
                                                  G4int          areacode[],
                                                  G4bool         isvalid[],
                                                  EValidate      validate = kValidateWithTol);
                                                  
   virtual G4int          DistanceToSurface(const G4ThreeVector &gp,
                                                  G4ThreeVector  gxx[],
                                                  G4double       distance[],
                                                  G4int          areacode[]);
 
   inline G4ThreeVector   ProjectAtPXPZ(const G4ThreeVector &p, G4bool isglobal = FALSE) const ;

 private:
   virtual G4double       DistanceToPlane(const G4ThreeVector &p,
                                          const G4ThreeVector &A,
                                          const G4ThreeVector &B,
                                          const G4ThreeVector &C,
                                          const G4ThreeVector &D,
                                          const G4int          parity,
                                                G4ThreeVector &xx,
                                                G4ThreeVector &n);

   virtual G4int          GetAreaCode(const G4ThreeVector &xx, 
                                            G4bool         withTol = TRUE);

   virtual void           SetCorners();
   virtual void           SetCorners(G4TwistedTubs *solid);
   virtual void           SetBoundaries();
   


private:

   G4double       fKappa;          // tan(TwistedAngle/2)/HalfLenZ;
  
};   


//========================================================
// inline function
//========================================================

inline
G4ThreeVector G4TwistedSurface::ProjectAtPXPZ(const G4ThreeVector &p, 
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
