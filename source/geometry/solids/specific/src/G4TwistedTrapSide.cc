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
// $Id: G4TwistedTrapSide.cc,v 1.3 2004-10-06 07:14:53 link Exp $
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
				     G4double      Halfzlen,
				     G4double      HalfSides[2],
				     G4double      AngleSide,
				     G4int         handedness)
  : G4VSurface(name)
{  

  fHandedness = handedness ;

  G4int ix = ( fHandedness > 0 ? 0 : 1 ) ; // ix = 0, iy = 1  if fHandedness > 0
  G4int iy = 1 - ix ;                      // ix = 1, iy = 1  if fHandedness < 0
	      
  fAxis[0]    = kYAxis; // in local coordinate system
  fAxis[1]    = kZAxis;
  fAxisMin[0] = -HalfSides[iy] ;  // Y Axis boundary
  fAxisMax[0] = HalfSides[iy] ;
  fAxisMin[1] = -Halfzlen ;      // Z Axis boundary
  fAxisMax[1] = Halfzlen ;
  
  fXHalfLength  = HalfSides[ix];  // b
  fYHalfLength  = HalfSides[iy];  // a
  fZHalfLength = Halfzlen ;  // L/2
  fPhiTwist = PhiTwist ;     // dphi
  fAngleSide = AngleSide ;  // 0,90,180,270 deg
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;

  // test
#ifdef G4SPECSDEBUG
  G4cout << "Side " << AngleSide << G4endl ;
  G4ThreeVector gp ( 100, 0, 0 ) ;
  G4ThreeVector p = ComputeLocalPoint(gp ) ;
  G4cout << "local point p = " << p << G4endl ;
#endif

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
   G4double phi = xx.z() / L  * fPhiTwist;
   G4double u   = ( fXHalfLength * cos(phi) - xx.x() ) / ( fXHalfLength * sin(phi) ) ;
   G4ThreeVector normal ( L * cos(phi) , L*sin(phi), fXHalfLength*fPhiTwist*u ) ;

#ifdef G4SPECSDEBUG
   G4cout  << "normal vector = " << normal << G4endl ;
   G4cout << "phi = " << phi << " , L = " << L << G4endl ;
#endif

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
      
  static const G4double ctol = 0.5 * kCarTolerance;

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
  
   // calculation of coefficients c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4 
  // with c4 = 1
 
  G4double L = 2*fZHalfLength ;
  G4double b = 2*fXHalfLength ;
  
  G4double phipzvx = fPhiTwist*p.z()*v.x()  ;
  G4double phipzvy = fPhiTwist*p.z()*v.y() ;
  G4double phipxvz = fPhiTwist*p.x()*v.z() ;
  G4double phipyvz = fPhiTwist*p.y()*v.z() ;

  G4double c[5],s[4] ;  

#ifdef G4SPECSDEBUG
  G4cout << " ---------------------------------- " << G4endl ;
  G4cout << " Point P = " << p << G4endl ;
  G4cout << " Direction v = " << v << G4endl ;
  G4cout << " b = " << b << G4endl ;
  G4cout << " L = " << L << G4endl ;
  G4cout << " twist = " << fPhiTwist << G4endl ;
  G4cout << " ---------------------------------- " << G4endl ;
#endif

  G4double ctmp = -4*L*v.y() + 6*b*fPhiTwist*v.z() ;
  c[0] = ( - 60*phipzvx - 30*b*fPhiTwist*v.z() + 60*phipxvz ) / ctmp  ;
  c[1] = ( 60*L*v.x() - 60*phipzvy + 60*phipyvz ) / ctmp  ;
  c[2] = ( 24*phipzvx + 60*L*v.y() - 3*b*fPhiTwist*v.z() - 24*phipxvz ) / ctmp ;
  c[3] = ( -24*L*v.x() + 4*phipzvy - 4*phipyvz ) / ctmp ;
  c[4] = 1 ;  // normal form of polynom (c4 = 1 )

#ifdef G4SPECSDEBUG
  G4double ct = ( - 60*phipzvx - 30*b*fPhiTwist*v.z() + 60*phipxvz ) ;
  G4cout << "c0 = " << ct << G4endl ;

  ct = 60*L*v.x()-60*phipzvy+60*phipyvz ;
  G4cout << "c1 = " << ct << G4endl ;

  ct = 24*phipzvx + 60*L*v.y() - 3*b*fPhiTwist*v.z() - 24*phipxvz ;
  G4cout << "c2 = " << ct << G4endl ;

  ct = -24*L*v.x() + 4*phipzvy - 4*phipyvz ;
  G4cout << "c3 = " << ct << G4endl ;

  ct = -4*L*v.y() + 6*b*fPhiTwist*v.z() ;
  G4cout << "c4 = " << ct << G4endl ;
#endif 


  // solve the polynom analytically
  G4ApproxPolySolver trapEq ;
  G4int num = trapEq.SolveBiQuadratic(c,s);

  // calculate phi and psi of the surface equation
  // and reduce the solution to the surface 

  G4double pi2 = 2*pi ;
  G4double tanPsi ;
  G4double PsiR = -1024  ;
  G4double PhiR = -1024 ;

  for (G4int i = 0 ; i<num ; i++ ) {

#ifdef G4SPECSDEBUG
    G4cout << "Solution " << i << " : " << s[i] << G4endl ;
#endif

    G4double stmp = fmod(s[i] , pi2) ;
    if ( s[i] < 0 ) { stmp -= 2*pi ; }
    G4double ztmp = L*s[i]/fPhiTwist ;
    if ( fabs(ztmp)<fZHalfLength+ctol ) {
      PhiR = stmp ;
      tanPsi = - 1 / ( b * fPhiTwist * v.z() ) 
	* ( 1./cos(PhiR) * 
	    ( 2*phipzvy - 2*phipyvz - 2*L*PhiR + b*fPhiTwist*v.z()*sin(PhiR) ) ) ;
      PsiR = atan(tanPsi) ;
#ifdef G4SPECSDEBUG
      G4cout << "solution " << i << " = " << PhiR << " , " << PsiR << " , " << ztmp << G4endl ;
#endif
    }
  }

  if ( PhiR < -1000 ) {
    //    G4cout << "TwistedTrapSide::DistanceToSurface problem with reconstruction" << G4endl ;
    distance[0] = kInfinity;
    gxx[0]      = ComputeGlobalPoint(xx[0]);
    isvalid[0] = false ;
    fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
				   areacode[0], isvalid[0],
				   0, validate, &gp, &gv);
    return 0;
  }

  // calculate the Point on the surface in cartesian coordinates
  // from the surface equations
  G4ThreeVector xxonsurface = SurfacePoint(PhiR,PsiR) ;

#ifdef G4SPECSDEBUG
  G4cout << "reconstructed phiR = " << PhiR << ", psiR = " << PsiR << G4endl ; 
#endif

  // the surfacenormal at that surface point
  G4ThreeVector surfacenormal = NormAng(PhiR,PsiR) ;

  // distance to that surfacepoint from particle position p along 
  // the direction v. The surface is approximated by a plane.
  distance[0] = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, xx[0]);

  // distance between reconstructed intersection point from the 
  // polynom-equation and the intersection point from distanceToPlaneWithV
  G4double deltaX = ( xx[0] - xxonsurface ).mag() ; 

  G4int maxint = 10 ;

  for ( G4int i = 1 ; i<maxint ; i++ ) {

    xxonsurface = SurfacePoint(PhiR,PsiR) ;
    surfacenormal = NormAng(PhiR,PsiR) ;
    distance[0] = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, xx[0]); // new XX[0]
    G4double deltaXtmp = ( xx[0] - xxonsurface ).mag() ; 

#ifdef G4SPECSDEBUG
    G4cout << "i = " << i << ", distance = " << distance[0] << ", " << deltaXtmp << G4endl ;
    G4cout << "X = " << xx[0] << G4endl ;
#endif
    if ( deltaX <= deltaXtmp && i> 1 ) { break ; } ;

    // the new point xx is accepted and phi/psi replaced
    GetPhiPsiAtX(xx[0], PhiR, PsiR) ;
    deltaX = deltaXtmp ;

    if ( deltaX <= 0.5*kCarTolerance ) { break ; }

    if ( i==maxint-1 ) {
      distance[0] = kInfinity;
      gxx[0]      = ComputeGlobalPoint(xx[0]);
      isvalid[0] = false ;
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
				     areacode[0], isvalid[0],
				     0, validate, &gp, &gv);
      return 0;
    }
      
  }

#ifdef G4SPECSDEBUG
  G4cout << "refined solution "  << PhiR << " , " << PsiR << " , " <<  G4endl ;
  G4cout << "distance = " << distance[0] << G4endl ;
  G4cout << "X = " << xx[0] << G4endl ;
#endif

  if (validate == kValidateWithTol) {
    areacode[0] = GetAreaCode(xx[0]);
    if (!IsOutside(areacode[0])) {
      if (distance[0] >= 0) isvalid[0] = true;
    }
  } else if (validate == kValidateWithoutTol) {
    areacode[0] = GetAreaCode(xx[0], false);
    if (IsInside(areacode[0])) {
      if (distance[0] >= 0) isvalid[0] = true;
    }
  } else { // kDontValidate
    G4Exception("G4TwistedTrapSide::DistanceToSurface()",
		"NotImplemented kDontValidate", FatalException,
		"Feature NOT implemented !");
  }

  gxx[0]      = ComputeGlobalPoint(xx[0]);

#ifdef G4SPECSDEBUG
  G4cout << "intersection Point found: " << gxx[0] << G4endl ;
  G4cout << "distance = " << distance[0] << G4endl ;
#endif

  fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
  return 1;

}



//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistedTrapSide::DistanceToSurface(const G4ThreeVector &gp,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[])
{  
  // to do

  static const G4double ctol = 0.5 * kCarTolerance;

  fCurStat.ResetfDone(kDontValidate, &gp);

   if (fCurStat.IsDone()) {
      G4int i;
      for (i=0; i<fCurStat.GetNXX(); i++) {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   } else {
      // initialize
      G4int i;
      for (i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector xx;  // intersection point
   G4ThreeVector xxonsurface ; // interpolated intersection point 

   // the surfacenormal at that surface point
   G4double phiR = 0  ; // 
   G4double psiR = 0 ;

   G4ThreeVector surfacenormal = NormAng(phiR,psiR) ;
   distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx);
   
   G4double deltaX = ( xx - xxonsurface ).mag() ; 
   
   G4int maxint = 10 ;

   for ( G4int i = 1 ; i<maxint ; i++ ) {

     xxonsurface = SurfacePoint(phiR,psiR) ;
     surfacenormal = NormAng(phiR,psiR) ;
     distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx); // new XX
     G4double deltaXtmp = ( xx - xxonsurface ).mag() ; 

#ifdef G4SPECSDEBUG
     G4cout << "i = " << i << ", distance = " << distance[0] << ", " << deltaXtmp << G4endl ;
     G4cout << "X = " << xx << G4endl ;
#endif

     if ( deltaX <= deltaXtmp && i> 1 ) { break ; } ;

    // the new point xx is accepted and phi/psi replaced
     GetPhiPsiAtX(xx, phiR, psiR) ;
     deltaX = deltaXtmp ;

     if ( deltaX <= ctol ) { break ; }

   }

   // check validity of solution ( valid phi,psi ) 

   G4double halfphi = 0.5*fPhiTwist ;
   G4double psiMax = atan(fYHalfLength/fXHalfLength) ;

   if (  phiR > halfphi ) phiR =  halfphi ;
   if ( phiR < -halfphi ) phiR = -halfphi ;
   if ( psiR > psiMax ) psiR = psiMax ;
   if ( psiR < -psiMax ) psiR = -psiMax ;

   xxonsurface = SurfacePoint(phiR,psiR) ;
   distance[0] = (  p - xx ).mag() ;
   if ( distance[0] <= ctol ) { distance[0] = 0 ; } 

   // end of validity 

#ifdef G4SPECSDEBUG
   G4cout << "refined solution "  << phiR << " , " << psiR << " , " <<  G4endl ;
   G4cout << "distance = " << distance[0] << G4endl ;
   G4cout << "X = " << xx << G4endl ;
#endif

   G4bool isvalid = true;
   gxx[0]      = ComputeGlobalPoint(xx);
   
#ifdef G4SPECSDEBUG
   G4cout << "intersection Point found: " << gxx[0] << G4endl ;
   G4cout << "distance = " << distance[0] << G4endl ;
#endif

   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
			    isvalid, 1, kDontValidate, &gp);
   return 1;
   

}


//=====================================================================
//* GetAreaCode -------------------------------------------------------

G4int G4TwistedTrapSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   static const G4double ctol = 0.5 * kCarTolerance;

   G4double phi = xx.z()/(2*fZHalfLength) * fPhiTwist ;

   G4double yprime = fXHalfLength * 
     ( fXHalfLength * cos(phi) - xx.x() ) / ( fXHalfLength * sin(phi) ) ; 

#ifdef G4SPECSDEBUG
   G4cout << "GetAreaCode: yprime = " << yprime << G4endl ;
#endif

   G4int areacode = sInside;
   
   if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
      G4int yaxis = 0;
      G4int zaxis = 1;
      
      if (withTol) {

	G4bool isoutside   = false;
	
	// test boundary of yaxis

         if (yprime < fAxisMin[yaxis] + ctol) {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary; 
            if (yprime <= fAxisMin[yaxis] - ctol) isoutside = true;

         } else if (yprime > fAxisMax[yaxis] - ctol) {
            areacode |= (sAxis0 & (sAxisY | sAxisMax)) | sBoundary;
            if (yprime >= fAxisMin[yaxis] + ctol)  isoutside = true;
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
            areacode |= (sAxis0 & sAxisY) | (sAxis1 & sAxisZ);
         }           
         
      } else {

         // boundary of y-axis

         if (yprime < fAxisMin[yaxis] ) {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary;
         } else if (yprime > fAxisMax[yaxis]) {
            areacode |= (sAxis0 & (sAxisY | sAxisMax)) | sBoundary;
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
            areacode |= (sAxis0 & sAxisY) | (sAxis1 & sAxisZ);
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

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
    
    G4double x, y, z;
    G4double bcos = fXHalfLength * cos(fPhiTwist/2) ;
    G4double bsin = fXHalfLength * sin(fPhiTwist/2) ;
    G4double acos = fYHalfLength * cos(fPhiTwist/2) ;
    G4double asin = fYHalfLength * sin(fPhiTwist/2) ;

    // corner of Axis0min and Axis1min
    x = bcos - asin ;
    y = -bsin - acos ;
    z = -fZHalfLength ;
    SetCorner(sCMin1Min, x, y, z);
      
    // corner of Axis0max and Axis1min
    x = bcos + asin ;
    y = -bsin + acos ;
    z = -fZHalfLength ;
    SetCorner(sCMax1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    x = bcos - asin ;
    y = bsin + acos ;
      z = fZHalfLength ;
    SetCorner(sCMax1Max, x, y, z);
      
    // corner of Axis0min and Axis1max
    x = bcos + asin ;
    y = bsin - acos ;
    z = fZHalfLength ;
    SetCorner(sCMin1Max, x, y, z);

  } else {

    G4Exception("G4TwistedTrapSide::SetCorners()",
		"NotImplemented", FatalException,
		"Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistedTrapSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   

  G4ThreeVector direction;
   
  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
      
    // sAxis0 & sAxisMin
    direction = GetCorner(sCMin1Max) - GetCorner(sCMin1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisY | sAxisMin), direction, 
		GetCorner(sCMin1Min), sAxisZ) ;
      
      // sAxis0 & sAxisMax
    direction = GetCorner(sCMax1Max) - GetCorner(sCMax1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisY | sAxisMax), direction, 
		GetCorner(sCMax1Min), sAxisZ);
    
    // sAxis1 & sAxisMin
    direction = GetCorner(sCMax1Min) - GetCorner(sCMin1Min);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
		GetCorner(sCMin1Min), sAxisY);
    
    // sAxis1 & sAxisMax
    direction = GetCorner(sCMax1Max) - GetCorner(sCMin1Max);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
		GetCorner(sCMin1Max), sAxisY);
    
  } else {
    
  G4Exception("G4TwistedTrapSide::SetCorners()",
	      "NotImplemented", FatalException,
	      "Feature NOT implemented !");
  }
  
}
