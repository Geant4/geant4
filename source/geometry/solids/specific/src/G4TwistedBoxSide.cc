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
// $Id: G4TwistedBoxSide.cc,v 1.1 2004-11-10 18:05:44 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedBoxSide.cc
//
// Author: 
//   Oliver Link
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistedBoxSide.hh"
//#include "G4PolynomialSolver.hh"
#include "G4ApproxPolySolver.hh" 

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistedBoxSide::G4TwistedBoxSide(const G4String     &name,
				   G4double      PhiTwist,
				   G4double      pDx,
				   G4double      pDy,
				   G4double      pDz,
				   G4double      AngleSide) 
  : G4VSurface(name)
{  

  // attention ! local coordinate system

  fAxis[0]    = kYAxis; // in local coordinate system
  fAxis[1]    = kZAxis;

  fDx = pDx ;
  fDy = pDy ;
  fDz = pDz ;

  fPhiTwist = PhiTwist ;     // dphi
  fAngleSide = AngleSide ;  // 0,90,180,270 deg
  
  fAxisMin[0] = -pDy ;  // Y Axis boundary
  fAxisMax[0] = pDy ;
  fAxisMin[1] = -pDz ;      // Z Axis boundary
  fAxisMax[1] = pDz ;
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;

}



//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedBoxSide::~G4TwistedBoxSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistedBoxSide::GetNormal(const G4ThreeVector &tmpxx, 
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

   G4double phi ;
   G4double u ;

   GetPhiUAtX(xx,phi,u) ;   // phi,u for point xx close to surface 

   G4ThreeVector normal =  NormAng(phi,u) ;  // the normal vector at phi,u

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

G4int G4TwistedBoxSide::DistanceToSurface(const G4ThreeVector &gp,
                                          const G4ThreeVector &gv,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[],
                                                G4bool         isvalid[],
                                                EValidate      validate)
{

  // implementation
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

  G4double uR = 0  ;
  G4double PhiR = -1024 ;
 
  G4double L = 2*fDz ;
  G4double a = 2*fDx ;

  // special case vz = 0

  if ( v.z() == 0. ) {         

    if ( fabs(p.z()) <= L ) {     // intersection possible in z
 
      PhiR = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 
      uR = (2*p.y()*v.x() - 2*p.x()*v.y() + a*v.y()*cos(PhiR) - a*v.x()*sin(PhiR) ) /
	(2*v.x()*cos(PhiR) + 2*v.y()*sin(PhiR)) ;
      //      G4cout  << "solution vz = 0 : "  << PhiR << " , " << uR << G4endl ;

    }

    else {                        // no intersection possible

      distance[0] = kInfinity;
      gxx[0]      = ComputeGlobalPoint(xx[0]);
      isvalid[0] = false ;
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
				     areacode[0], isvalid[0],
				     0, validate, &gp, &gv);
      return 0;


    }  // end fabs(p.z() <= L 

  } // end v.z() == 0
  

  // general solution for non-zero vz

  else {

    G4double phipzvx = fPhiTwist*p.z()*v.x()  ;
    G4double phipzvy = fPhiTwist*p.z()*v.y() ;

    G4double phivz   = fPhiTwist*v.z();
    G4double phipxvz = phivz * p.x() ;
    G4double phipyvz = phivz * p.y() ;

    G4double c[5],s[4] ;  

    // calculation of coefficients c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4 
    // with c4 = 1

#ifdef G4SPECSDEBUG
    G4cout << " ---------------------------------- " << G4endl ;
    G4cout << " Point P = " << p << G4endl ;
    G4cout << " Direction v = " << v << G4endl ;
    G4cout << " a = " << a << G4endl ;
    G4cout << " L = " << L << G4endl ;
    G4cout << " twist = " << fPhiTwist << G4endl ;
    G4cout << " ---------------------------------- " << G4endl ;
#endif
    
    G4double ctmp = -4*L*v.y() + 6*a*phivz ;
    c[0] = ( - 60*phipzvx - 30*a*phivz + 60*phipxvz ) / ctmp  ;
    c[1] = ( 60*L*v.x() - 60*phipzvy + 60*phipyvz ) / ctmp  ;
    c[2] = ( 24*phipzvx + 60*L*v.y() - 3*a*phivz - 24*phipxvz ) / ctmp ;
    c[3] = ( -24*L*v.x() + 4*phipzvy - 4*phipyvz ) / ctmp ;
    c[4] = 1 ;  // normal form of polynom (c4 = 1 )
    
  // solve the polynom analytically
    G4ApproxPolySolver trapEq ;
    G4int num = trapEq.SolveBiQuadratic(c,s);
  
  // calculate phi and psi of the surface equation
  // and reduce the solution to the surface 

    G4double pi2 = 2*pi ;

    for (G4int i = 0 ; i<num ; i++ ) {

#ifdef G4SPECSDEBUG
      G4cout << "Solution " << i << " : " << s[i] << G4endl ;
#endif

      G4double stmp = fmod(s[i] , pi2) ;
      if ( s[i] < 0 && stmp > 0 ) { stmp -= 2*pi ; }
      G4double ztmp = L*s[i]/fPhiTwist ;
      if ( fabs(ztmp)<fDz+ctol ) {
	PhiR = stmp ;
	uR = (2*L*PhiR*v.y() - 2*phipzvy + 2*phipyvz - a*phivz * sin(PhiR))/(2.*phivz*cos(PhiR)) ;
	
#ifdef G4SPECSDEBUG
	G4cout << "solution " << i << " = " << PhiR << " , " << uR  << G4endl ;
#endif
      }
    }

    if ( PhiR < -1000 ) {
    //    G4cout << "TwistedBoxSide::DistanceToSurface problem with reconstruction" << G4endl ;
      distance[0] = kInfinity;
      gxx[0]      = ComputeGlobalPoint(xx[0]);
      isvalid[0] = false ;
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
				     areacode[0], isvalid[0],
				     0, validate, &gp, &gv);
      return 0;

    }

  }    // end general case



  // calculate the Point on the surface in cartesian coordinates
  // from the surface equations
  G4ThreeVector xxonsurface = SurfacePoint(PhiR,uR) ;

#ifdef G4SPECSDEBUG
  G4cout << "reconstructed phiR = " << PhiR << ", uR = " << uR << G4endl ; 
#endif

  // the surfacenormal at that surface point
  G4ThreeVector surfacenormal = NormAng(PhiR,uR) ;

  // distance to that surfacepoint from particle position p along 
  // the direction v. The surface is approximated by a plane.
  distance[0] = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, xx[0]);

  // distance between reconstructed intersection point from the 
  // polynom-equation and the intersection point from distanceToPlaneWithV
  G4double deltaX = ( xx[0] - xxonsurface ).mag() ; 

  G4int maxint = 10 ;

  for ( G4int i = 1 ; i<maxint ; i++ ) {

    xxonsurface = SurfacePoint(PhiR,uR) ;
    surfacenormal = NormAng(PhiR,uR) ;
    distance[0] = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, xx[0]); // new XX[0]
    G4double deltaXtmp = ( xx[0] - xxonsurface ).mag() ; 

#ifdef G4SPECSDEBUG
    G4cout << "i = " << i << ", distance = " << distance[0] << ", " << deltaXtmp << G4endl ;
    G4cout << "X = " << xx[0] << G4endl ;
#endif
    if ( deltaX <= deltaXtmp && i> 1 ) { break ; } ;

    // the new point xx is accepted and phi/psi replaced
    GetPhiUAtX(xx[0], PhiR, uR) ;
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
  G4cout << "refined solution "  << PhiR << " , " << uR << " , " <<  G4endl ;
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
    G4Exception("G4TwistedBoxSide::DistanceToSurface()",
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

G4int G4TwistedBoxSide::DistanceToSurface(const G4ThreeVector &gp,
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
   G4double uR = 0 ;

   G4ThreeVector surfacenormal = NormAng(phiR,uR) ;
   distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx);
   
   G4double deltaX = ( xx - xxonsurface ).mag() ; 
   
   G4int maxint = 10 ;

   for ( G4int i = 1 ; i<maxint ; i++ ) {

     xxonsurface = SurfacePoint(phiR,uR) ;
     surfacenormal = NormAng(phiR,uR) ;
     distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx); // new XX
     G4double deltaXtmp = ( xx - xxonsurface ).mag() ; 

#ifdef G4SPECSDEBUG
     G4cout << "i = " << i << ", distance = " << distance[0] << ", " << deltaXtmp << G4endl ;
     G4cout << "X = " << xx << G4endl ;
#endif

     if ( deltaX <= deltaXtmp && i> 1 ) { break ; } ;

    // the new point xx is accepted and phi/psi replaced
     GetPhiUAtX(xx, phiR, uR) ;
     deltaX = deltaXtmp ;

     if ( deltaX <= ctol ) { break ; }

   }

   // check validity of solution ( valid phi,psi ) 

   G4double halfphi = 0.5*fPhiTwist ;
   G4double uMax = fDy ;

   if (  phiR > halfphi ) phiR =  halfphi ;
   if ( phiR < -halfphi ) phiR = -halfphi ;
   if ( uR > uMax ) uR = uMax ;
   if ( uR < -uMax ) uR = -uMax ;

   xxonsurface = SurfacePoint(phiR,uR) ;
   distance[0] = (  p - xx ).mag() ;
   if ( distance[0] <= ctol ) { distance[0] = 0 ; } 

   // end of validity 

#ifdef G4SPECSDEBUG
   G4cout << "refined solution "  << phiR << " , " << uR << " , " <<  G4endl ;
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

G4int G4TwistedBoxSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   static const G4double ctol = 0.5 * kCarTolerance;

   G4double phi ;
   G4double yprime ;
   GetPhiUAtX(xx, phi,yprime ) ;

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
      G4Exception("G4TwistedBoxSide::GetAreaCode()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistedBoxSide::SetCorners()
{

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
    
    G4double x, y, z;
    G4double acos = fDx * cos(fPhiTwist/2) ;
    G4double asin = fDx * sin(fPhiTwist/2) ;
    G4double bcos = fDy * cos(fPhiTwist/2) ;
    G4double bsin = fDy * sin(fPhiTwist/2) ;

    // corner of Axis0min and Axis1min
    x = acos - bsin ;
    y = -asin - bcos ;
    z = -fDz ;
    SetCorner(sC0Min1Min, x, y, z);
      
    // corner of Axis0max and Axis1min
    x = acos + bsin ;
    y = -asin + bcos ;
    z = -fDz ;
    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    x = acos - bsin ;
    y = asin + bcos ;
      z = fDz ;
    SetCorner(sC0Max1Max, x, y, z);
      
    // corner of Axis0min and Axis1max
    x = acos + bsin ;
    y = asin - bcos ;
    z = fDz ;
    SetCorner(sC0Min1Max, x, y, z);

  } else {

    G4Exception("G4TwistedBoxSide::SetCorners()",
		"NotImplemented", FatalException,
		"Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistedBoxSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   

  G4ThreeVector direction;
   
  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
      
    // sAxis0 & sAxisMin
    direction = GetCorner(sC0Min1Max) - GetCorner(sC0Min1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisY | sAxisMin), direction, 
		GetCorner(sC0Min1Min), sAxisZ) ;
      
      // sAxis0 & sAxisMax
    direction = GetCorner(sC0Max1Max) - GetCorner(sC0Max1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisY | sAxisMax), direction, 
		GetCorner(sC0Max1Min), sAxisZ);
    
    // sAxis1 & sAxisMin
    direction = GetCorner(sC0Max1Min) - GetCorner(sC0Min1Min);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
		GetCorner(sC0Min1Min), sAxisY);
    
    // sAxis1 & sAxisMax
    direction = GetCorner(sC0Max1Max) - GetCorner(sC0Min1Max);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
		GetCorner(sC0Min1Max), sAxisY);
    
  } else {
    
  G4Exception("G4TwistedBoxSide::SetCorners()",
	      "NotImplemented", FatalException,
	      "Feature NOT implemented !");
  }
  
}



void G4TwistedBoxSide::GetPhiUAtX( G4ThreeVector p, G4double &phi, G4double &u) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction ( both for a fixed z = pz)
  
  phi = p.z()/(2*fDz)*fPhiTwist ;
  G4ThreeVector  X0 ( fDx * cos(phi), fDx * sin(phi) , p.z() ) ;  // basis with u=0
  G4ThreeVector dvec  ( - sin(phi), cos(phi) ,0. ) ;   // direction vector
  G4ThreeVector xx ;                                   // the intersection point on the line

  DistanceToLine(p ,X0 , dvec , xx) ;
  
  u = ( xx - X0 ).mag() ;  // X0 is choosen such that u = 0

}


G4ThreeVector G4TwistedBoxSide::ProjectPoint(const G4ThreeVector &p, 
                                                    G4bool isglobal) 
{
  // Get Rho at p.z() on Hyperbolic Surface.
  G4ThreeVector tmpp;
  if (isglobal) {
     tmpp = fRot.inverse()*p - fTrans;
  } else {
     tmpp = p;
  }

  G4double phi ;
  G4double u ;

  GetPhiUAtX( tmpp, phi, u ) ;  // calculate (phi, u) for a point p close the surface
  
  G4ThreeVector xx = SurfacePoint(phi,u) ;  // transform back to cartesian coordinates

  if (isglobal) {
     return (fRot * xx + fTrans);
  } else {
     return xx;
  }
}


