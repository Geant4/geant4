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
// $Id: G4TwistedTrapSide.cc,v 1.1 2004-07-29 15:10:49 link Exp $
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
				     G4int         handedness)
  : G4VSurface(name)
{  
   fHandedness = handedness;   // +z = +ve, -z = -ve
   fAxis[0]    = kXAxis; // in local coordinate system
   fAxis[1]    = kZAxis;
   fAxisMin[0] = InnerRadius;  // Inner-hype radius at z=0
   fAxisMax[0] = OuterRadius;  // Outer-hype radius at z=0
   fAxisMin[1] = EndZ[0];
   fAxisMax[1] = EndZ[1];


   fKappa = Kappa;
   fSideX  = 2*HalfSides[0];
   fSideY  = 2*HalfSides[1];
   fZHalfLength = halfzlen ;
   fPhiTwist = PhiTwist ;

   G4cout << "zlen/2 = " << fZHalfLength << G4endl ;
   G4cout << "PhiTwist = " << fPhiTwist << G4endl ;
   

   fRot.rotateZ( fHandedness > 0
                 ? -0.5*DPhi
                 :  0.5*DPhi );
   fTrans.set(0, 0, 0);
   fIsValidNorm = false;
   
   SetCorners( EndInnerRadius, EndOuterRadius, EndPhi, EndZ) ;
   SetBoundaries();
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
   normal = normal/normal.mag() ;

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

  G4ThreeVector p = gp;
  G4ThreeVector v = gv;
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
   fCurStat.ResetfDone(kDontValidate, &gp);

   if (fCurStat.IsDone()) {
      for (G4int i=0; i<fCurStat.GetNXX(); i++) {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   } else {
      // initialize
      for (G4int i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   static const G4double halftol = 0.5 * kCarTolerance; 

   G4ThreeVector  p       = ComputeLocalPoint(gp);
   G4ThreeVector  xx;
   G4int          parity  = (fKappa >= 0 ? 1 : -1);
 
   // 
   // special case! 
   // If p is on surface, or
   // p is on z-axis, 
   // return here immediatery.
   //
   
   G4ThreeVector  lastgxx[2];
   G4double       distfromlast[2];
   for (G4int i=0; i<2; i++) {
      lastgxx[i] = fCurStatWithV.GetXX(i);
      distfromlast[i] = (gp - lastgxx[i]).mag();
   } 

   if  ((gp - lastgxx[0]).mag() < halftol
     || (gp - lastgxx[1]).mag() < halftol) { 
      // last winner, or last poststep point is on the surface.
      xx = p;
      distance[0] = 0;
      gxx[0] = gp;

      G4bool isvalid = true;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
      return 1;
   }
          
   if (p.getRho() == 0) { 
      // p is on z-axis. Namely, p is on twisted surface (invalid area).
      // We must return here, however, returning distance to x-minimum
      // boundary is better than return 0-distance.
      //
      G4bool isvalid = true;
      if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
         distance[0] = DistanceToBoundary(sAxis0 & sAxisMin, xx, p);
         areacode[0] = sInside;
      } else {
         distance[0] = 0;
         xx.set(0., 0., 0.);
      }
      gxx[0] = ComputeGlobalPoint(xx);
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 0, kDontValidate, &gp);
      return 1;
   } 

   // 
   // special case end
   //

   // set corner points of quadrangle try area ...

   G4ThreeVector A;  // foot of normal from p to boundary of sAxis0 & sAxisMin
   G4ThreeVector C;  // foot of normal from p to boundary of sAxis0 & sAxisMax
   G4ThreeVector B;       // point on boundary sAxis0 & sAxisMax at z = A.z()
   G4ThreeVector D;       // point on boundary sAxis0 & sAxisMin at z = C.z()
   G4double      distToA; // distance from p to A
   G4double      distToC; // distance from p to C 

   distToA = DistanceToBoundary(sAxis0 & sAxisMin, A, p);
   distToC = DistanceToBoundary(sAxis0 & sAxisMax, C, p);
   
   // is p.z between a.z and c.z?
   // p.z must be bracketed a.z and c.z.
   if (A.z() > C.z()) {
      if (p.z() > A.z()) {
         A = GetBoundaryAtPZ(sAxis0 & sAxisMin, p);
      } else if (p.z() < C.z()) {
         C = GetBoundaryAtPZ(sAxis0 & sAxisMax, p);
      }
   } else {
      if (p.z() > C.z()) {
         C = GetBoundaryAtPZ(sAxis0 & sAxisMax, p);
      } else if (p.z() < A.z()) {
         A = GetBoundaryAtPZ(sAxis0 & sAxisMin, p);
      }
   }

   G4ThreeVector  d[2];     // direction vectors of boundary
   G4ThreeVector  x0[2];    // foot of normal from line to p 
   G4int          btype[2]; // boundary type

   for (G4int i=0; i<2; i++) {
      if (i == 0) {
         GetBoundaryParameters((sAxis0 & sAxisMax), d[i], x0[i], btype[i]);
         B = x0[i] + ((A.z() - x0[i].z()) / d[i].z()) * d[i]; 
         // x0 + t*d , d is direction unit vector.
      } else {
         GetBoundaryParameters((sAxis0 & sAxisMin), d[i], x0[i], btype[i]);
         D = x0[i] + ((C.z() - x0[i].z()) / d[i].z()) * d[i]; 
      }
   }

   // In order to set correct diagonal, swap A and D, C and B if needed.  
   G4ThreeVector pt(p.x(), p.y(), 0.);
   G4double      rc = fabs(p.x());
   G4ThreeVector surfacevector(rc, rc * fKappa * p.z(), 0.); 
   G4int         pside = AmIOnLeftSide(pt, surfacevector); 
   G4double      test  = (A.z() - C.z()) * parity * pside;  

   if (test == 0) {
      if (pside == 0) {
         // p is on surface.
         xx = p;
         distance[0] = 0;
         gxx[0] = gp;

         G4bool isvalid = true;
         fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);
         return 1;
      } else {
         // A.z = C.z(). return distance to line.
         d[0] = C - A;
         distance[0] = DistanceToLine(p, A, d[0], xx);
         areacode[0] = sInside;
         gxx[0] = ComputeGlobalPoint(xx);
         G4bool isvalid = true;
         fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);
         return 1;
      } 

   } else if (test < 0) {

      // wrong diagonal. vector AC is crossing the surface!  
      // swap A and D, C and B
      G4ThreeVector tmp;
      tmp = A;
      A = D;
      D = tmp;
      tmp = C;
      C = B;
      B = tmp; 

   } else {
      // correct diagonal. nothing to do.  
   }

   // Now, we chose correct diaglnal.
   // First try. divide quadrangle into double triangle by diagonal and 
   // calculate distance to both surfaces.

   G4ThreeVector xxacb;   // foot of normal from plane ACB to p
   G4ThreeVector nacb;    // normal of plane ACD
   G4ThreeVector xxcad;   // foot of normal from plane CAD to p
   G4ThreeVector ncad;    // normal of plane CAD
   G4ThreeVector AB(A.x(), A.y(), 0);
   G4ThreeVector DC(C.x(), C.y(), 0);

   G4double distToACB = G4VSurface::DistanceToPlane(p, A, C-A, AB, xxacb, nacb) * parity;
   G4double distToCAD = G4VSurface::DistanceToPlane(p, C, C-A, DC, xxcad, ncad) * parity;

   // if calculated distance = 0, return  

   if (fabs(distToACB) <= halftol || fabs(distToCAD) <= halftol) {
      xx = (fabs(distToACB) < fabs(distToCAD) ? xxacb : xxcad); 
      areacode[0] = sInside;
      gxx[0] = ComputeGlobalPoint(xx);
      distance[0] = 0;
      G4bool isvalid = true;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0] , areacode[0],
                                isvalid, 1, kDontValidate, &gp);
      return 1;
   }
   
   if (distToACB * distToCAD > 0 && distToACB < 0) {
      // both distToACB and distToCAD are negative.
      // divide quadrangle into double triangle by diagonal
      G4ThreeVector normal;
      distance[0] = DistanceToPlane(p, A, B, C, D, parity, xx, normal);
   } else {
      if (distToACB * distToCAD > 0) {
         // both distToACB and distToCAD are positive.
         // Take smaller one.
         if (distToACB <= distToCAD) {
            distance[0] = distToACB;
            xx   = xxacb;
         } else {
            distance[0] = distToCAD;
            xx   = xxcad;
         }
      } else {
         // distToACB * distToCAD is negative.
         // take positive one
         if (distToACB > 0) {
            distance[0] = distToACB;
            xx   = xxacb;
         } else {
            distance[0] = distToCAD;
            xx   = xxcad;
         }
      }
      
   }
   areacode[0] = sInside;
   gxx[0]      = ComputeGlobalPoint(xx);
   G4bool isvalid = true;
   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
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
   static const G4double halftol = 0.5 * kCarTolerance;
   
   G4ThreeVector M = 0.5*(A + B);
   G4ThreeVector N = 0.5*(C + D);
   G4ThreeVector xxanm;  // foot of normal from p to plane ANM
   G4ThreeVector nanm;   // normal of plane ANM
   G4ThreeVector xxcmn;  // foot of normal from p to plane CMN
   G4ThreeVector ncmn;   // normal of plane CMN

   G4double distToanm = G4VSurface::DistanceToPlane(p, A, (N - A), (M - A), xxanm, nanm) * parity;
   G4double distTocmn = G4VSurface::DistanceToPlane(p, C, (M - C), (N - C), xxcmn, ncmn) * parity;

   // if p is behind of both surfaces, abort.
   if (distToanm * distTocmn > 0 && distToanm < 0) {
     G4Exception("G4TwistedTrapSide::DistanceToPlane()",
                 "InvalidCondition", FatalException,
                 "Point p is behind the surfaces.");
   }

   // if p is on surface, return 0.
   if (fabs(distToanm) <= halftol) {
      xx = xxanm;
      n  = nanm * parity;
      return 0;
   } else if (fabs(distTocmn) <= halftol) {
      xx = xxcmn;
      n  = ncmn * parity;
      return 0;
   }
   
   if (distToanm <= distTocmn) {
      if (distToanm > 0) {
         // both distanses are positive. take smaller one.
         xx = xxanm;
         n  = nanm * parity;
         return distToanm;
      } else {
         // take -ve distance and call the function recursively.
         return DistanceToPlane(p, A, M, N, D, parity, xx, n);
      }
   } else {
      if (distTocmn > 0) {
         // both distanses are positive. take smaller one.
         xx = xxcmn;
         n  = ncmn * parity;
         return distTocmn;
      } else {
         // take -ve distance and call the function recursively.
         return DistanceToPlane(p, C, N, M, B, parity, xx, n);
      }
   }
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
//* SetCorners( arglist ) -------------------------------------------------

void G4TwistedTrapSide::SetCorners(
                                  G4double      endInnerRad[2],
                                  G4double      endOuterRad[2],
                                  G4double      endPhi[2],
                                  G4double      endZ[2])
{
   // Set Corner points in local coodinate.   

   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
   
      G4int zmin = 0 ;  // at -ve z
      G4int zmax = 1 ;  // at +ve z

      G4double x, y, z;
      
      // corner of Axis0min and Axis1min
      x = endInnerRad[zmin]*cos(endPhi[zmin]);
      y = endInnerRad[zmin]*sin(endPhi[zmin]);
      z = endZ[zmin];
      SetCorner(sCMin1Min, x, y, z);
      
      // corner of Axis0max and Axis1min
      x = endOuterRad[zmin]*cos(endPhi[zmin]);
      y = endOuterRad[zmin]*sin(endPhi[zmin]);
      z = endZ[zmin];
      SetCorner(sCMax1Min, x, y, z);
      
      // corner of Axis0max and Axis1max
      x = endOuterRad[zmax]*cos(endPhi[zmax]);
      y = endOuterRad[zmax]*sin(endPhi[zmax]);
      z = endZ[zmax];
      SetCorner(sCMax1Max, x, y, z);
      
      // corner of Axis0min and Axis1max
      x = endInnerRad[zmax]*cos(endPhi[zmax]);
      y = endInnerRad[zmax]*sin(endPhi[zmax]);
      z = endZ[zmax];
      SetCorner(sCMin1Max, x, y, z);

   } else {
      G4cerr << "ERROR - G4FlatSurface::SetCorners()" << G4endl
             << "        fAxis[0] = " << fAxis[0] << G4endl
             << "        fAxis[1] = " << fAxis[1] << G4endl;
      G4Exception("G4TwistedTrapSide::SetCorners()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
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
   G4ThreeVector direction;
   
   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
      
      // sAxis0 & sAxisMin
      direction = GetCorner(sCMin1Max) - GetCorner(sCMin1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisX | sAxisMin), direction, 
                  GetCorner(sCMin1Min), sAxisZ) ;
      
      // sAxis0 & sAxisMax
      direction = GetCorner(sCMax1Max) - GetCorner(sCMax1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisX | sAxisMax), direction, 
                  GetCorner(sCMax1Min), sAxisZ);
                  
      // sAxis1 & sAxisMin
      direction = GetCorner(sCMax1Min) - GetCorner(sCMin1Min);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
                  GetCorner(sCMin1Min), sAxisX);
                  
      // sAxis1 & sAxisMax
      direction = GetCorner(sCMax1Max) - GetCorner(sCMin1Max);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
                  GetCorner(sCMin1Max), sAxisX);
                  
   } else {
      G4cerr << "ERROR - G4FlatSurface::SetBoundaries()" << G4endl
             << "        fAxis[0] = " << fAxis[0] << G4endl
             << "        fAxis[1] = " << fAxis[1] << G4endl;
      G4Exception("G4TwistedTrapSide::SetCorners()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
}

