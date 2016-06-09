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
// $Id: G4TwistedTrapSide.cc,v 1.16 2005/02/14 13:55:52 link Exp $
// GEANT4 tag $Name: geant4-07-00-patch-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrapSide.cc
//
// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// History:
//    14.2.05 Changed Polynom Solver to G4JTPolynomialSolver
//
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistedTrapSide.hh"
#include "G4JTPolynomialSolver.hh"

//=====================================================================
//* constructors ------------------------------------------------------


G4TwistedTrapSide::G4TwistedTrapSide(const G4String     &name,
                                     G4double      PhiTwist,
                                     G4double      pDx1,
                                     G4double      pDx2,
                                     G4double      pDy,
                                     G4double      pDz,
                                     G4double      AngleSide)
  : G4VSurface(name)
{  

              
  fAxis[0]    = kYAxis; // in local coordinate system
  fAxis[1]    = kZAxis;
  fAxisMin[0] = -pDy ;  // Y Axis boundary
  fAxisMax[0] = pDy ;
  fAxisMin[1] = -pDz ;      // Z Axis boundary
  fAxisMax[1] = pDz ;
  
  fDx1  = pDx1 ;
  fDx2  = pDx2 ;
  fDy   = pDy ;
  fDz   = pDz ;

  fPhiTwist = PhiTwist ;     // dphi
  fAngleSide = AngleSide ;  // 0,90,180,270 deg
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;

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

   G4double phi ;
   G4double u ;

   GetPhiUAtX(xx,phi,u) ;   // phi,u for point xx close to surface 

   G4ThreeVector normal =  NormAng(phi,u) ;  // the normal vector at phi,u

#ifdef G4SPECSDEBUG
   G4cout  << "normal vector = " << normal << G4endl ;
   G4cout << "phi = " << phi << " , u = " << u << G4endl ;
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
{

  static const G4double ctol = 0.5 * kCarTolerance;
  static const G4double pihalf = pi/2 ;

  G4bool IsParallel = false ;
  G4bool IsConverged =  false ;

  G4int nxx = 0 ;  // number of physical solutions

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
    for (i=0; i<G4VSURFACENXX ; i++) {
      distance[i] = kInfinity;
      areacode[i] = sOutside;
      isvalid[i]  = false;
      gxx[i].set(kInfinity, kInfinity, kInfinity);
    }
  }

  G4ThreeVector p = ComputeLocalPoint(gp);
  G4ThreeVector v = ComputeLocalDirection(gv);
  
#ifdef G4SPECSDEBUG
  G4cout << "Local point p = " << p << G4endl ;
  G4cout << "Local direction v = " << v << G4endl ; 
#endif

  G4double phi,u ;  // parameters

  // temporary variables

  G4double      tmpdist = kInfinity ;
  G4ThreeVector tmpxx;
  G4int         tmpareacode = sOutside ;
  G4bool        tmpisvalid  = false ;

  std::vector<Intersection> xbuf ;
  Intersection xbuftmp ;
  
  // prepare some variables for the intersection finder

  G4double L = 2*fDz ;
  G4double a = 2*fDx2 ;
  G4double d = 2*fDx1 ;
  G4double b = 2*fDy ;

  G4double amd = a-d ;
  G4double apd = a+d ;

  G4double phipzvx = fPhiTwist*p.z()*v.x() ;
  G4double phipzvy = fPhiTwist*p.z()*v.y() ;

  G4double phipyvz = fPhiTwist*p.y()*v.z() ;
  G4double phipxvz = fPhiTwist*p.x()*v.z() ;

  G4double phivz   = fPhiTwist*v.z();
  G4double phiyz =   fPhiTwist * ( v.z()*p.y() - v.y()*p.z() )   ;

  // special case vz = 0

  if ( v.z() == 0. ) {         

    if ( std::fabs(p.z()) <= L ) {     // intersection possible in z
      
      phi = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 
      u =  (b*(4*p.y()*v.x() - 4*p.x()*v.y() + apd*v.y()*std::cos(phi) 
	       - apd*v.x()*std::sin(phi)))/
        (2.*((2*b*v.x() - amd *v.y())*std::cos(phi) 
	     + ( amd * v.x() + 2*b*v.y())*std::sin(phi)));

      xbuftmp.phi = phi ;
      xbuftmp.u = u ;
      xbuftmp.areacode = sOutside ;
      xbuftmp.distance = kInfinity ;
      xbuftmp.isvalid = false ;

      xbuf.push_back(xbuftmp) ;  // store it to xbuf

    }

    else {                        // no intersection possible

      distance[0] = kInfinity;
      gxx[0].set(kInfinity,kInfinity,kInfinity);
      isvalid[0] = false ;
      areacode[0] = sOutside ;
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
                                     areacode[0], isvalid[0],
                                     0, validate, &gp, &gv);
      
      return 0;


    }  // end std::fabs(p.z() <= L 

  } // end v.z() == 0
  

  // general solution for non-zero vz

  else {

    G4double c[8],sr[7],si[7] ;  


#ifdef G4SPECSDEBUG
    G4cout << " ---------------------------------- " << G4endl ;
    G4cout << " Point P = " << p << G4endl ;
    G4cout << " Direction v = " << v << G4endl ;
    G4cout << " a = " << a << G4endl ;
    G4cout << " d = " << d << G4endl ;
    G4cout << " b = " << b << G4endl ;
    G4cout << " L = " << L << G4endl ;
    G4cout << " twist = " << fPhiTwist << G4endl ;
    G4cout << " ---------------------------------- " << G4endl ;
#endif
    
    c[7] =  3600*fPhiTwist*(b*(4*p.z()*v.x() + (a + d - 4*p.x())*v.z()) - 2*amd*(p.z()*v.y() - p.y()*v.z())) ;
    c[6] = -7200*(-(amd*(phipzvx + L*v.y() - phipxvz)) + 2*b*(L*v.x() - phipzvy + phipyvz)) ;
    c[5] = 300*(-4*b*(5*phipzvx + 12*L*v.y() - 5*phipxvz) + 
		d*(24*L*v.x() - 10*phipzvy + b*phivz + 10*phipyvz) + 
		a*(-24*L*v.x() + fPhiTwist*(10*p.z()*v.y() + b*v.z() - 10*p.y()*v.z()))) ;
    c[4] = 600*(-(amd*(phipzvx + 5*L*v.y() - phipxvz)) +  2*b*(5*L*v.x() - phipzvy + phipyvz)) ;
    c[3] = 6*(100*amd*L*v.x() + 4*b*phipzvx + 200*b*L*v.y() - 2*amd*phipzvy + 
	      a*b*fPhiTwist*v.z() + b*d*fPhiTwist*v.z() - 4*b*phipxvz + 2*amd*phipyvz) ;
    c[2] = 4*(-(amd*(7*phipzvx - 3*L*v.y() - 7*phipxvz)) - 2*b*(3*L*v.x() + 7*phipzvy - 7*phipyvz)) ;
    c[1] = 28*amd*L*v.x() + 18*b*phipzvx + 56*b*L*v.y() - 9*amd*phipzvy - 18*b*phipxvz + 9*amd*phipyvz ;
    c[0] = -9*L*(2*b*v.x() -amd*v.y()) ;

#ifdef G4SPECSDEBUG
    G4cout << "coef = " << c[0] << " " 
	   <<  c[1] << " "  
	   <<  c[2] << " "  
	   <<  c[3] << " "  
	   <<  c[4] << " "  
	   <<  c[5] << " "  
	   <<  c[6] << " "  
	   <<  c[7] << G4endl ;
#endif    

    G4JTPolynomialSolver trapEq ;

    G4int num = trapEq.FindRoots(c,7,sr,si);
  

    for (G4int i = 0 ; i<num ; i++ ) {  // loop over all mathematical solutions
      if ( si[i]==0.0 ) {  // only real solutions
#ifdef G4SPECSDEBUG
	G4cout << "Solution " << i << " : " << sr[i] << G4endl ;
#endif
	phi = std::fmod(sr[i] , pihalf)  ;

	u   = b * ( 4*phiyz + 4*L*v.y()*phi - apd*phivz*std::sin(phi) ) 
	  / (2.*phivz*(2*b*std::cos(phi) + amd*std::sin(phi) )) ;
	
	xbuftmp.phi = phi ;
	xbuftmp.u = u ;
	xbuftmp.areacode = sOutside ;
	xbuftmp.distance = kInfinity ;
	xbuftmp.isvalid = false ;
	
	xbuf.push_back(xbuftmp) ;  // store it to xbuf
      
#ifdef G4SPECSDEBUG
	G4cout << "solution " << i << " = " << phi << " , " << u  << G4endl ;
#endif

      }  // end if real solution
    }  // end loop i
    
  }    // end general case


  nxx = xbuf.size() ;  // save the number of  solutions

  G4ThreeVector xxonsurface  ;       // point on surface
  G4ThreeVector surfacenormal  ;     // normal vector  
  G4double deltaX  ;                 // distance between intersection point and point on surface
  G4double theta  ;                  // angle between track and surfacenormal
  G4double factor ;                  // a scaling factor
  G4int maxint = 30 ;                // number of iterations


  for ( size_t k = 0 ; k<xbuf.size() ; k++ ) {

#ifdef G4SPECSDEBUG
    G4cout << "Solution " << k << " : " 
	   << "reconstructed phiR = " << xbuf[k].phi
	   << ", uR = " << xbuf[k].u << G4endl ; 
#endif
    
    phi = xbuf[k].phi ;  // get the stored values for phi and u
    u = xbuf[k].u ;

    IsConverged = false ;   // no convergence at the beginning
    
    for ( G4int i = 1 ; i<maxint ; i++ ) {
      
      xxonsurface = SurfacePoint(phi,u) ;
      surfacenormal = NormAng(phi,u) ;
      tmpdist = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, tmpxx); 
      deltaX = ( tmpxx - xxonsurface ).mag() ; 
      theta = std::fabs(std::acos(v*surfacenormal) - pihalf) ;
      if ( theta < 0.001 ) { 
	factor = 50 ;
	IsParallel = true ;
      }
      else {
	factor = 1 ;
      }

#ifdef G4SPECSDEBUG
      G4cout << "Step i = " << i << ", distance = " << tmpdist << ", " << deltaX << G4endl ;
      G4cout << "X = " << tmpxx << G4endl ;
#endif
      
      GetPhiUAtX(tmpxx, phi, u) ; // the new point xx is accepted and phi/u replaced
      
#ifdef G4SPECSDEBUG
      G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
      if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
      
    }  // end iterative loop (i)
    

#ifdef G4SPECSDEBUG
    G4cout << "refined solution "  << phi << " , " << u  <<  G4endl ;
    G4cout << "distance = " << tmpdist << G4endl ;
    G4cout << "local X = " << tmpxx << G4endl ;
#endif
    
    tmpisvalid = false ;  // init 

    if ( IsConverged ) {
      
      if (validate == kValidateWithTol) {
	tmpareacode = GetAreaCode(tmpxx);
	if (!IsOutside(tmpareacode)) {
	  if (tmpdist >= 0) tmpisvalid = true;
	}
      } else if (validate == kValidateWithoutTol) {
	tmpareacode = GetAreaCode(tmpxx, false);
	if (IsInside(tmpareacode)) {
	  if (tmpdist >= 0) tmpisvalid = true;
	}
      } else { // kDontValidate
	G4Exception("G4TwistedTrapSide::DistanceToSurface()",
		    "NotImplemented kDontValidate", FatalException,
		    "Feature NOT implemented !");
      }

    } 
    else {
      tmpdist = kInfinity;     // no convergence after 10 steps 
      tmpisvalid = false ;     // solution is not vaild
    }  


    // store the found values 
    xbuf[k].xx = tmpxx ;
    xbuf[k].distance = tmpdist ;
    xbuf[k].areacode = tmpareacode ;
    xbuf[k].isvalid = tmpisvalid ;


  }  // end loop over physical solutions (variable k)


  std::sort(xbuf.begin() , xbuf.end(), DistanceSort ) ;  // sorting

#ifdef G4SPECSDEBUG
  G4cout << G4endl << "list xbuf after sorting : " << G4endl ;
  G4cout << G4endl << G4endl ;
#endif


  // erase identical intersection (within kCarTolerance) 
  xbuf.erase( std::unique(xbuf.begin(), xbuf.end() , EqualIntersection ) , xbuf.end() ) ;


  // add guesses

  G4int nxxtmp = xbuf.size() ;

  if ( nxxtmp<2 || IsParallel  ) {

    // positive end
#ifdef G4SPECSDEBUG
    G4cout << "add guess at +z/2 .. " << G4endl ;
#endif

    phi = fPhiTwist/2 ;
    u   = b * ( 4*phiyz + 4*L*v.y()*phi - apd*phivz*std::sin(phi) ) 
      / (2.*phivz*(2*b*std::cos(phi) + amd*std::sin(phi) )) ;
    
    xbuftmp.phi = phi ;
    xbuftmp.u = u ;
    xbuftmp.areacode = sOutside ;
    xbuftmp.distance = kInfinity ;
    xbuftmp.isvalid = false ;
    
    xbuf.push_back(xbuftmp) ;  // store it to xbuf


#ifdef G4SPECSDEBUG
    G4cout << "add guess at -z/2 .. " << G4endl ;
#endif

    phi = -fPhiTwist/2 ;
    u   = b * ( 4*phiyz + 4*L*v.y()*phi - apd*phivz*std::sin(phi) ) 
      / (2.*phivz*(2*b*std::cos(phi) + amd*std::sin(phi) )) ;
    
    xbuftmp.phi = phi ;
    xbuftmp.u = u ;
    xbuftmp.areacode = sOutside ;
    xbuftmp.distance = kInfinity ;
    xbuftmp.isvalid = false ;
    
    xbuf.push_back(xbuftmp) ;  // store it to xbuf

    for ( size_t k = nxxtmp ; k<xbuf.size() ; k++ ) {

#ifdef G4SPECSDEBUG
      G4cout << "Solution " << k << " : " 
	     << "reconstructed phiR = " << xbuf[k].phi
	     << ", uR = " << xbuf[k].u << G4endl ; 
#endif
      
      phi = xbuf[k].phi ;  // get the stored values for phi and u
      u   = xbuf[k].u ;

      IsConverged = false ;   // no convergence at the beginning
      
      for ( G4int i = 1 ; i<maxint ; i++ ) {
	
	xxonsurface = SurfacePoint(phi,u) ;
	surfacenormal = NormAng(phi,u) ;
	tmpdist = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, tmpxx); 
	deltaX = ( tmpxx - xxonsurface ).mag() ; 
	theta = std::fabs(std::acos(v*surfacenormal) - pihalf) ;
	if ( theta < 0.001 ) { 
	  factor = 50 ;    
	}
	else {
	  factor = 1 ;
	}
	
#ifdef G4SPECSDEBUG
	G4cout << "Step i = " << i << ", distance = " << tmpdist << ", " << deltaX << G4endl ;
	G4cout << "X = " << tmpxx << G4endl ;
#endif

	GetPhiUAtX(tmpxx, phi, u) ; // the new point xx is accepted and phi/u replaced
      
#ifdef G4SPECSDEBUG
	G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
	if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
      
      }  // end iterative loop (i)
    

#ifdef G4SPECSDEBUG
      G4cout << "refined solution "  << phi << " , " << u  <<  G4endl ;
      G4cout << "distance = " << tmpdist << G4endl ;
      G4cout << "local X = " << tmpxx << G4endl ;
#endif

      tmpisvalid = false ;  // init 

      if ( IsConverged ) {

	if (validate == kValidateWithTol) {
	  tmpareacode = GetAreaCode(tmpxx);
	  if (!IsOutside(tmpareacode)) {
	    if (tmpdist >= 0) tmpisvalid = true;
	  }
	} else if (validate == kValidateWithoutTol) {
	  tmpareacode = GetAreaCode(tmpxx, false);
	  if (IsInside(tmpareacode)) {
	    if (tmpdist >= 0) tmpisvalid = true;
	  }
	} else { // kDontValidate
	  G4Exception("G4TwistedBoxSide::DistanceToSurface()",
		      "NotImplemented kDontValidate", FatalException,
		      "Feature NOT implemented !");
	}
	
      } 
      else {
	tmpdist = kInfinity;     // no convergence after 10 steps 
	tmpisvalid = false ;     // solution is not vaild
      }  
	
	
      // store the found values 
      xbuf[k].xx = tmpxx ;
      xbuf[k].distance = tmpdist ;
      xbuf[k].areacode = tmpareacode ;
      xbuf[k].isvalid = tmpisvalid ;


    }  // end loop over physical solutions 


  }  // end less than 2 solutions


  // sort again
  std::sort(xbuf.begin() , xbuf.end(), DistanceSort ) ;  // sorting

  // erase identical intersection (within kCarTolerance) 
  xbuf.erase( std::unique(xbuf.begin(), xbuf.end() , EqualIntersection ) , xbuf.end() ) ;

#ifdef G4SPECSDEBUG
  G4cout << G4endl << "list xbuf after sorting : " << G4endl ;
  G4cout << G4endl << G4endl ;
#endif

  nxx = xbuf.size() ;   // determine number of solutions again.

  for ( size_t i = 0 ; i<xbuf.size() ; i++ ) {
    
    distance[i] = xbuf[i].distance;
    gxx[i]      = ComputeGlobalPoint(xbuf[i].xx);
    areacode[i] = xbuf[i].areacode ;
    isvalid[i]  = xbuf[i].isvalid ;
    
    fCurStatWithV.SetCurrentStatus(i, gxx[i], distance[i], areacode[i],
				     isvalid[i], nxx, validate, &gp, &gv);

#ifdef G4SPECSDEBUG
    G4cout << "element Nr. " << i 
	   << ", local Intersection = " << xbuf[i].xx 
	   << ", distance = " << xbuf[i].distance 
	   << ", u = " << xbuf[i].u 
	   << ", phi = " << xbuf[i].phi 
	   << ", isvalid = " << xbuf[i].isvalid 
	   << G4endl ;
#endif

  }  // end for( i ) loop

    
#ifdef G4SPECSDEBUG
  G4cout << "G4TwistedTrapSide finished " << G4endl ;
  G4cout << nxx << " possible physical solutions found" << G4endl ;
  for ( G4int k= 0 ; k< nxx ; k++ ) {
    G4cout << "global intersection Point found: " << gxx[k] << G4endl ;
    G4cout << "distance = " << distance[k] << G4endl ;
    G4cout << "isvalid = " << isvalid[k] << G4endl ;
  }
#endif

  return nxx ;
    
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
      for (i=0; i<G4VSURFACENXX; i++) {
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

   G4ThreeVector surfacenormal ; 
   G4double deltaX ;
   
   G4int maxint = 20 ;

   for ( G4int i = 1 ; i<maxint ; i++ ) {

     xxonsurface = SurfacePoint(phiR,uR) ;
     surfacenormal = NormAng(phiR,uR) ;
     distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx); // new XX
     deltaX = ( xx - xxonsurface ).mag() ; 

#ifdef G4SPECSDEBUG
     G4cout << "i = " << i << ", distance = " << distance[0] << ", " << deltaX << G4endl ;
     G4cout << "X = " << xx << G4endl ;
#endif

     // the new point xx is accepted and phi/psi replaced
     GetPhiUAtX(xx, phiR, uR) ;
     
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

G4int G4TwistedTrapSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   static const G4double ctol = 0.5 * kCarTolerance;

   G4double phi ;
   G4double yprime ;
   GetPhiUAtX(xx, phi,yprime ) ;

#ifdef G4SPECSDEBUG
   G4cout << "GetAreaCode: phi = " << phi << G4endl ;
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
    G4double a1cos = fDx1 * std::cos(fPhiTwist/2) ;
    G4double a1sin = fDx1 * std::sin(fPhiTwist/2) ;
    G4double a2cos = fDx2 * std::cos(fPhiTwist/2) ;
    G4double a2sin = fDx2 * std::sin(fPhiTwist/2) ;
    G4double bcos = fDy * std::cos(fPhiTwist/2) ;
    G4double bsin = fDy * std::sin(fPhiTwist/2) ;

    // corner of Axis0min and Axis1min
    x = a1cos - bsin ;
    y = -a1sin - bcos ;
    z = -fDz ;
    SetCorner(sC0Min1Min, x, y, z);
      
    // corner of Axis0max and Axis1min
    x = a2cos + bsin ;
    y = -a2sin + bcos ;
    z = -fDz ;
    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    x = a2cos - bsin ;
    y = a2sin + bcos ;
      z = fDz ;
    SetCorner(sC0Max1Max, x, y, z);
      
    // corner of Axis0min and Axis1max
    x = a1cos + bsin ;
    y = a1sin - bcos ;
    z = fDz ;
    SetCorner(sC0Min1Max, x, y, z);

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
    
  G4Exception("G4TwistedTrapSide::SetCorners()",
              "NotImplemented", FatalException,
              "Feature NOT implemented !");
  }
  
}



void G4TwistedTrapSide::GetPhiUAtX( G4ThreeVector p, G4double &phi, G4double &u) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction ( both for a fixed z = pz)
  
#if 0
  phi = p.z()/(2*fDz)*fPhiTwist ;
  G4ThreeVector  X0 ( Xcoef(0.) * std::cos(phi), Xcoef(0.) * std::sin(phi) , p.z() ) ;  // basis with u=0
  G4ThreeVector dvec  ( - (fDx1-fDx2)/(2*fDy) * std::cos(phi) - std::sin(phi), 
                        std::cos(phi) - (fDx1-fDx2)/(2*fDy)*std::sin(phi) ,
                        0. ) ;   // direction vector

  G4ThreeVector xx ;                                   // the intersection point on the line

  DistanceToLine(p ,X0 , dvec , xx) ;
  
  u = ( xx - X0 ).mag() / dvec.mag() ;  // X0 is choosen such that u = 0

#endif

  // phi is given by the z coordinate of p

  phi = p.z()/(2*fDz)*fPhiTwist ;

  G4double cphi = std::cos(phi) ;
  G4double sphi = std::sin(phi) ;
  G4double c0 = (fDx2-fDx1)/(2*fDy) ;

  // this formula is the analytical form of the procedure used above.
  // this is not faster, but removes some spourious events.
  // See Mathematica notebook: dirTrap.m 

  u = ( (fDx1*fDx1-fDx2*fDx2)/(4*fDy) + c0*cphi*p.x() + c0*sphi*p.y() + 
	( p.y() * cphi - p.x() * sphi ) ) / 
    ( ( c0*cphi - sphi)*(c0*cphi - sphi)  + ( cphi + c0*sphi )*(cphi + c0*sphi ));

}


G4ThreeVector G4TwistedTrapSide::ProjectPoint(const G4ThreeVector &p, 
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


