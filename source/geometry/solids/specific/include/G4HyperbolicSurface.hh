// $Id: G4HyperbolicSurface.hh,v 1.1 2003-11-12 15:43:00 link Exp $
#ifndef __G4HYPERBOLICSURFACE__
#define __G4HYPERBOLICSURFACE__
//*************************************************************************
//* --------------------
//* G4HyperbolicSurface
//* --------------------
//* (Description)
//* 	Class for hyperbolic boundary surface for J4Solid.
//*     
//* (Update Record)
//*	2002/08/01  K.Hoshina	Original version.
//*************************************************************************

#include "G4VSurface.hh"

class G4TwistedTubs;

class G4HyperbolicSurface : public G4VSurface
{
public:
                       
   G4HyperbolicSurface(const G4String         &name,
                       const G4RotationMatrix &rot,       // rot of 0.5*(phi-width of a segment)
                       const G4ThreeVector    &tlate,
                       const G4int             handedness,// right hand = 1, left hand = -1
                       const G4double          kappa,     // tan(TwistAngle/2)/fZHalfLen
                       const G4double          tanstereo, // tan(stereo angle)
                       const G4double          r0,        // radius at z = 0
                       const EAxis             axis0 = kPhi,
                       const EAxis             axis1 = kZAxis,
                             G4double         axis0min = -kInfinity,
                             G4double         axis1min = -kInfinity,
                             G4double         axis0max = kInfinity,
                             G4double         axis1max = kInfinity); 
                             
   G4HyperbolicSurface(const G4String      &name,
                             G4TwistedTubs *solid,
                             G4int          handedness);  // right hand = 1, left hand = -1

  
   virtual ~G4HyperbolicSurface();

   virtual G4int           DistanceToSurface(const G4ThreeVector &gp,
                                             const G4ThreeVector &gv,
                                                   G4ThreeVector  gxx[],
                                                   G4double       distance[],
                                                   G4int          areacode[],
                                                   G4bool         isvalid[],
                                                   EValidate      validate = kValidateWithTol);
                                                   
   virtual G4int           DistanceToSurface(const G4ThreeVector &gp,
                                                   G4ThreeVector  gxx[],
                                                   G4double       distance[],
                                                   G4int          areacode[]);
 
   virtual G4ThreeVector   GetNormal(const G4ThreeVector &xx, G4bool isGlobal = FALSE) ;
   virtual EInside         Inside(const G4ThreeVector &gp) ;
   
   inline virtual G4double GetRhoAtPZ(const G4ThreeVector &p, G4bool isglobal = FALSE) const ;
   
private:
   virtual G4int           GetAreaCode(const G4ThreeVector &xx, 
                                             G4bool withTol = TRUE);
   virtual G4int           GetAreaCodeInPhi(const G4ThreeVector &xx, 
                                                  G4bool withTol = TRUE);

   virtual void            SetCorners();
   virtual void            SetCorners(G4TwistedTubs *solid);
   virtual void            SetBoundaries();

private:
   
   G4double          fKappa;	      // tan(TwistedAngle/2)/HalfLenZ;
   G4double          fTanStereo;    // tan(StereoAngle)
   G4double          fTan2Stereo;   // tan(StereoAngle)**2
   G4double          fR0;           // radius at z = 0
   G4double          fR02;          // radius**2 at z = 0
   class             Insidetype{
                       public:
                         G4ThreeVector gp;
                         EInside       inside;
                     };
   Insidetype        fInside;

};


//========================================================
// inline function
//========================================================

inline
G4double G4HyperbolicSurface::GetRhoAtPZ(const G4ThreeVector &p, G4bool isglobal) const 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }
  return sqrt(fR02 + tmpp.z() * tmpp.z() * fTan2Stereo); 
}

#endif

