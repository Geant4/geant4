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
// $Id: G4TwistedTrapSide.cc,v 1.2 2004-08-27 13:36:59 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrapSide.cc
//
// Author: 
//   Oliver Link
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistedTrapSide.hh"
//#include "G4PolynomialSolver.hh"
#include "G4ApproxPolySolver.hh" 

//=====================================================================
//* constructors ------------------------------------------------------


G4TwistedTrapSide::G4TwistedTrapSide(const G4String     &name,
				     G4double      PhiTwist,
				     G4double      halfzlen,
				     G4double      HalfSides[2],
				     G4int         handedness)
  : G4VSurface(name)
{  
   fHandedness = handedness;   // +z = +ve, -z = -ve
   fAxis[0]    = kXAxis; // in local coordinate system
   fAxis[1]    = kZAxis;

   fSideX  = 2*HalfSides[0];
   fSideY  = 2*HalfSides[1];
   fZHalfLength = halfzlen ;
   fPhiTwist = PhiTwist ;

   fRot.rotateZ( fHandedness ) ; 

   fTrans.set(0, 0, 0);
   fIsValidNorm = false;
   

}



//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedTrapSide::~G4TwistedTrapSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistedTrapSide::GetNormal(const G4ThreeVector &tmpxx, 
                                                G4bool isGlobal) 
{
   // GetNormal returns a normal vector at a surface (or very close
   // to surface) point at tmpxx.
   // If isGlobal=true, it returns the normal in global coordinate.
   //
   G4ThreeVector xx;
   if (isGlobal) {
      xx = ComputeLocalPoint(tmpxx);
      if ((xx - fCurrentNormal.p).mag() < 0.5 * kCarTolerance) {
         return ComputeGlobalDirection(fCurrentNormal.normal);
      }
   } else {
      xx = tmpxx;
      if (xx == fCurrentNormal.p) {
         return fCurrentNormal.normal;
      }
   }

   G4double L   = 2*fZHalfLength ;   
   G4double phi = tmpxx.z() / L  * fPhiTwist;
   G4double u   = ( fSideY/2 * cos(phi) - tmpxx.x() ) / ( fSideY/2 * sin(phi) ) ;
   G4ThreeVector normal = ( L * cos(phi) , L*sin(phi), fSideY/2*fPhiTwist*u ) ;

   //    normal = normal/normal.mag() ;

   if (isGlobal) {
      fCurrentNormal.normal = ComputeGlobalDirection(normal.unit());
   } else {
      fCurrentNormal.normal = normal.unit();
   }
   return fCurrentNormal.normal;
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistedTrapSide::DistanceToSurface(const G4ThreeVector &gp,
                                          const G4ThreeVector &gv,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[],
                                                G4bool         isvalid[],
                                                EValidate      validate)
{  // implementation
   // Coordinate system:
   //
      
  //  G4cout << "particle position  = " << gp << G4endl ;
  // G4cout << "particle direction = " << gv << G4endl ;

  fCurStatWithV.ResetfDone(validate, &gp, &gv);

  if (fCurStatWithV.IsDone()) {
    G4int i;
    for (i=0; i<fCurStatWithV.GetNXX(); i++) {
      gxx[i] = fCurStatWithV.GetXX(i);
      distance[i] = fCurStatWithV.GetDistance(i);
      areacode[i] = fCurStatWithV.GetAreacode(i);
      isvalid[i]  = fCurStatWithV.IsValid(i);
    }
    return fCurStatWithV.GetNXX();
  } else {
   
   // initialize
    G4int i;
    for (i=0; i<2; i++) {
      distance[i] = kInfinity;
      areacode[i] = sOutside;
      isvalid[i]  = false;
      gxx[i].set(kInfinity, kInfinity, kInfinity);
    }
  }

  G4ThreeVector p = ComputeLocalPoint(gp);
  G4ThreeVector v = ComputeLocalDirection(gv);

  G4ThreeVector xx[2]; 
  
  // 
  // special case! 
  // p is origin or
  //
  
  G4double absvz = fabs(v.z());

  if ((absvz < DBL_MIN) && (fabs(p.x() * v.y() - p.y() * v.x()) < DBL_MIN)) {
      // no intersection
    
    isvalid[0] = false;
    fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
			      isvalid[0], 0, validate, &gp, &gv);
    return 0;
  } 
   
  // 
  // special case end
  //
  
  // calculation of coefficients c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4 
  
  G4double L = 2*fZHalfLength ;
  G4double b = 2*fSideY ;
  
  G4double phipzvx = fPhiTwist*p.z()*v.x()  ;
  G4double phipzvy = fPhiTwist*p.z()*v.y() ;
  G4double phipxvz = fPhiTwist*p.x()*v.z() ;
  G4double phipyvz = fPhiTwist*p.y()*v.z() ;

  G4double c[5],s[4] ;  

  G4double ctmp = -4*L*v.y() + 6*b*fPhiTwist*v.z() ;
  c[0] = ( - 60*phipzvx - 30*b*fPhiTwist*v.z() + 60*phipxvz ) / ctmp  ;
  c[1] = ( 60*L*v.x()-60*phipxvz+60*phipyvz ) / ctmp  ;
  c[2] = ( 24*phipzvx + 60*L*v.y() - 3*b*fPhiTwist*v.z() - 24*phipxvz ) / ctmp ;
  c[3] = ( -24*L*v.x() + 4*phipzvy - 4*phipyvz ) / ctmp ;
  c[4] = 1 ;  // normal form of polynom (c4 = 1 )

  // solve the polynom analytically
  G4ApproxPolySolver trapEq ;
  G4int num = trapEq.SolveBiQuadratic(c,s);


  // calculate phi and psi of the surface equation
  G4double pi2 = 2*pi ;
  G4double tanPsi ;
  G4double PsiR = -1024 ;
  G4double PhiR = -1024 ;

  for (G4int i = 0 ; i<num ; i++ ) {
    G4double stmp = fmod(s[i] , pi2) ;
    if ( s[i] < 0 ) { stmp -= 2*pi ; }
    //    G4cout << "solution " << i << " = " << s[i] 
    //	   << " -> stmp = " << stmp << G4endl ;
    if ( fabs(stmp)<fPhiTwist/2.) { // good solution 
      PhiR = stmp ;
      tanPsi = - 1 / ( b * fPhiTwist * v.z() ) 
	* ( 1./cos(PhiR) * 
	    ( 2*phipzvy - 2*phipyvz - 2*L*PhiR + b*fPhiTwist*v.z()*sin(PhiR) ) ) ;
      PsiR = atan(tanPsi) ;
    }
  }

  if ( PhiR > -1000 ) {
    //    G4cout << "reconstructed phiR = " << PhiR << ", psiR = " << PsiR << G4endl ; 

 
  // calculate the Point on the surface in cartesian coordinates
  // from the surface equations

    G4ThreeVector xxonsurface = SurfacePoint(PhiR,PsiR) ;

    // the surfacenormal at that surface point
    
    G4ThreeVector surfacenormal = NormAng(PhiR,PsiR) ;

    // distance to that surfacepoint from particle position p along 
    // the direction v. The surface is approximated by a plane.
    distance[0] = DistanceToPlaneWithV(p, v, xxonsurface,
				       surfacenormal, xx[0]);

    // G4cout << "distance = " << distance[0] << G4endl ;

// distance between reconstructed intersection point from the 
// polynom-equation and the intersection point from distanceToPlaneWithV
    G4double deltaX = ( xx[0] - xxonsurface ).mag() ; 

    G4int maxint = 10 ;

    for ( G4int i = 1 ; i<maxint ; i++ ) {

      xxonsurface = SurfacePoint(PhiR,PsiR) ;
      surfacenormal = NormAng(PhiR,PsiR) ;
      distance[0] = DistanceToPlaneWithV(p, v, xxonsurface,
				       surfacenormal, xx[0]);
      G4double deltaXtmp = ( xx[0] - xxonsurface ).mag() ; 

      //      G4cout << " i = " << i << ", distance = " << distance[0] << ", " << deltaXtmp << G4endl ;
      // G4cout << "X = " << xx[0] << G4endl ;
      if ( deltaX <= deltaXtmp && i> 1 ) { break ; } ;

      // the new point xx is accepted and phi/psi replaced
      GetPhiPsiAtX(xx[0], PhiR, PsiR) ;
      deltaX = deltaXtmp ;

      if ( deltaX <= 0.5*kCarTolerance ) { break ; }

      if ( i==maxint-1 ) {
	distance[0] = kInfinity;
	gxx[0]      = ComputeGlobalPoint(xx[0]);
	fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
				       areacode[0], isvalid[0],
				       0, validate, &gp, &gv);
	return 0;
      }

    }
    

    for ( G4int i = 0 ; i < 1 ; i++ ) {
      if (validate == kValidateWithTol) {
	areacode[i] = GetAreaCode(xx[i]);
	if (!IsOutside(areacode[i])) {
	  if (distance[i] >= 0) isvalid[i] = true;
	  continue;
	}
      } else if (validate == kValidateWithoutTol) {
	areacode[i] = GetAreaCode(xx[i], false);
	if (IsInside(areacode[i])) {
	  if (distance[i] >= 0) isvalid[i] = true;
	  continue;
	}
      } else { // kDontValidate
	// we must choose x(rho,z) = rho(z=0) > 0
	if (xx[i].x() > 0) {
	  areacode[i] = sInside;
	  if (distance[i] >= 0) isvalid[i] = true;
	  continue;
	} else {
	  distance[i] = kInfinity;
	  continue;
	}                     
      }
    }
  
    gxx[0]      = ComputeGlobalPoint(xx[0]);

    fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
    return 1;

  } else {
    return 0 ;
  }

}



//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistedTrapSide::DistanceToSurface(const G4ThreeVector &gp,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[])
{  
  // to do
   
   return 1;
}

//=====================================================================
//* DistanceToPlane ---------------------------------------------------

G4double G4TwistedTrapSide::DistanceToPlane(const G4ThreeVector &p,
                                           const G4ThreeVector &A,
                                           const G4ThreeVector &B,
                                           const G4ThreeVector &C,
                                           const G4ThreeVector &D,
                                           const G4int          parity,
                                                 G4ThreeVector &xx,
                                                 G4ThreeVector &n)
{

  // to do
  return 0.0 ;
}

//=====================================================================
//* GetAreaCode -------------------------------------------------------

G4int G4TwistedTrapSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   static const G4double ctol = 0.5 * kCarTolerance;
   G4int areacode = sInside;
   
   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
      G4int xaxis = 0;
      G4int zaxis = 1;
      
      if (withTol) {

         G4bool isoutside   = false;

         // test boundary of xaxis

         if (xx.x() < fAxisMin[xaxis] + ctol) {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary; 
            if (xx.x() <= fAxisMin[xaxis] - ctol) isoutside = true;

         } else if (xx.x() > fAxisMax[xaxis] - ctol) {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
            if (xx.x() >= fAxisMin[xaxis] + ctol)  isoutside = true;
         }

         // test boundary of z-axis

         if (xx.z() < fAxisMin[zaxis] + ctol) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin)); 

            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary;
            if (xx.z() <= fAxisMin[zaxis] - ctol) isoutside = true;

         } else if (xx.z() > fAxisMax[zaxis] - ctol) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax));

            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
            if (xx.z() >= fAxisMax[zaxis] + ctol) isoutside = true;
         }

         // if isoutside = true, clear inside bit.             
         // if not on boundary, add axis information.             
         
         if (isoutside) {
            G4int tmpareacode = areacode & (~sInside);
            areacode = tmpareacode;
         } else if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
         
      } else {

         // boundary of x-axis

         if (xx.x() < fAxisMin[xaxis] ) {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary;
         } else if (xx.x() > fAxisMax[xaxis]) {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
         }
         
         // boundary of z-axis

         if (xx.z() < fAxisMin[zaxis]) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
           
         } else if (xx.z() > fAxisMax[zaxis]) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax)) ;
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
         }

         if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
      }
      return areacode;
   } else {
      G4Exception("G4TwistedTrapSide::GetAreaCode()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistedTrapSide::SetCorners()
{
   G4Exception("G4TwistedTrapSide::SetCorners()",
               "NotImplemented", FatalException,
               "Method NOT implemented !");
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistedTrapSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   
  G4Exception("G4TwistedTrapSide::SetCorners()",
	      "NotImplemented", FatalException,
	      "Feature NOT implemented !");
}

