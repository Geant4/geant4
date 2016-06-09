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
// $Id: G4TwistedTrapAlphaSide.cc,v 
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrapAlphaSide.cc
//
// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistedTrapAlphaSide.hh"
#include "G4JTPolynomialSolver.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistedTrapAlphaSide::G4TwistedTrapAlphaSide(const G4String     &name,
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
			   G4double      pAlph,      // tilt angle at +pDz
                           G4double      AngleSide    // parity
					       ) : G4VSurface(name)
{  
  
                 
  fAxis[0]    = kYAxis; // in local coordinate system
  fAxis[1]    = kZAxis;
  fAxisMin[0] = -kInfinity ;  // Y Axis boundary
  fAxisMax[0] = kInfinity ;   //   depends on z !!
  fAxisMin[1] = -pDz ;      // Z Axis boundary
  fAxisMax[1] = pDz ;
  
  fDx1  = pDx1 ;
  fDx2  = pDx2 ;
  fDx3  = pDx3 ;
  fDx4  = pDx4 ;

  fDy1   = pDy1 ;
  fDy2   = pDy2 ;

  fDz   = pDz ;

  fAlph = pAlph  ;
  fTAlph = std::tan(fAlph) ;

  fTheta = pTheta ;
  fPhi   = pPhi ;

  // precalculate frequently used parameters
  fDx4plus2  = fDx4 + fDx2 ;
  fDx4minus2 = fDx4 - fDx2 ;
  fDx3plus1  = fDx3 + fDx1 ; 
  fDx3minus1 = fDx3 - fDx1 ;
  fDy2plus1  = fDy2 + fDy1 ;
  fDy2minus1 = fDy2 - fDy1 ;

  fa1md1 = 2*fDx2 - 2*fDx1  ; 
  fa2md2 = 2*fDx4 - 2*fDx3 ;

  fPhiTwist = PhiTwist ;     // dphi
  fAngleSide = AngleSide ;  // 0,90,180,270 deg

  fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;  // dx in surface equation
  fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;  // dy in surface equation
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;

}


//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedTrapAlphaSide::~G4TwistedTrapAlphaSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistedTrapAlphaSide::GetNormal(const G4ThreeVector &tmpxx, 
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

G4int G4TwistedTrapAlphaSide::DistanceToSurface(const G4ThreeVector &gp,
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

  //  G4double phipzvx = fPhiTwist*p.z()*v.x() ;
  // G4double phipzvy = fPhiTwist*p.z()*v.y() ;

  //  G4double phipxvz = fPhiTwist*p.x()*v.z() ;
  // G4double phipyvz = fPhiTwist*p.y()*v.z() ;

  // G4double phivz   = fPhiTwist*v.z();
  G4double phix    = fPhiTwist * ( v.x()*p.z() - v.z()*p.x() )   ;
  G4double phiy    = fPhiTwist * ( v.y()*p.z() - v.z()*p.y() )   ;


  // special case vz = 0

  if ( v.z() == 0. ) {         

    if ( std::fabs(p.z()) <= L ) {     // intersection possible in z
      
      phi = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 
      u = (2*fDy1*(8*fPhiTwist*(p.y()*v.x() - p.x()*v.y()) + 
          (2*(fDx4plus2 + fDx3plus1)*fPhiTwist*v.y() - 2*(4*fdeltaY*v.x() + (-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*v.y())*phi)*std::cos(phi) - 
          (2*(fDx4plus2 + fDx3plus1)*fPhiTwist*v.x() - 2*((-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*v.x() - 4*fdeltaY*v.y())*phi)*std::sin(phi)))/
      (4.*fPhiTwist*(std::sin(phi)*(fa1md1*v.x() + 4*fDy1*v.y() - 4*fDy1*v.x()*fTAlph) + std::cos(phi)*(4*fDy1*v.x() - fa1md1 *v.y()  + 4*fDy1*v.y()*fTAlph)))  ;

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

    c[7] = -1800*(-4*fa1md1*phiy + 16*fDy1*(phix + phiy*fTAlph) +  4*(fDx4plus2 + fDx3plus1)*fDy1*fPhiTwist*v.z()) ;
    c[6] = -3600*(8*fDy1*(phiy - L*v.x()) + 2*(phix - fdeltaY*v.z() + L*v.y())*(fa1md1 - 4*fDy1*fTAlph) - 
		  2*fDy1*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*v.z()) ;
    c[5] = 150*(fa1md1*(-20*phiy + 48*L*v.x()) + 2*fDy1*(40*phix + 96*L*v.y() + 40*phiy*fTAlph - 96*L*v.x()*fTAlph - 
							 2*(fDx4plus2 + fDx3plus1)*fPhiTwist*v.z())) ;
    c[4] = -300*(-8*fDy1*(phiy - 5*L*v.x()) - 2*(phix + fdeltaY*v.z() + 5*L*v.y())*(fa1md1 - 4*fDy1*fTAlph) - 
		 2*fDy1*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*v.z()) ;
    c[3] = -3*(-4*fa1md1*(phiy - 50*L*v.x()) + 16*fDy1*(phix + 50*L*v.y() + phiy*fTAlph - 50*L*v.x()*fTAlph) + 
	       2*(fDx4plus2 + fDx3plus1)*2*fDy1*fPhiTwist*v.z()) ;
    c[2] = 16*fDy1*(7*phiy + 3*L*v.x()) + 4*(7*phix + 3*fdeltaY*v.z() 
					     - 3*L*v.y())*(fa1md1 - 4*fDy1*fTAlph) + 12*fDy1*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*v.z() ;
    c[1] = -4*fDy1*(9*fPhiTwist*p.z()*v.x() - 9*fPhiTwist*p.x()*v.z() + 
		    28*L*v.y()) + (9*phiy - 28*L*v.x())*(fa1md1 - 4*fDy1*fTAlph)  ;
    c[0] = 9*(-(fa1md1*L*v.y()) + 4*fDy1*(L*v.x() + L*v.y()*fTAlph)) ;

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

	u   = (16*fDy1*(-phiy + L*v.y()*phi - fdeltaY*v.z()*phi*std::cos(phi)) - 
     2*fDy1*v.z()*( 2*(fDx4plus2 + fDx3plus1)*fPhiTwist -2*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*phi)*std::sin(phi))/
   (4.*fPhiTwist*v.z()*(4*fDy1*(std::cos(phi) - std::sin(phi)*fTAlph) + fa1md1*std::sin(phi)))  ;
	
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
	G4Exception("G4TwistedTrapAlphaSide::DistanceToSurface()",
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
    u   =  (16*fDy1*(-phiy + L*v.y()*phi - fdeltaY*v.z()*phi*std::cos(phi)) - 
     2*fDy1*v.z()*( 2*(fDx4plus2 + fDx3plus1)*fPhiTwist -2*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*phi)*std::sin(phi))/
   (4.*fPhiTwist*v.z()*(4*fDy1*(std::cos(phi) - std::sin(phi)*fTAlph) + fa1md1*std::sin(phi))) ;
    
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
    u   = (16*fDy1*(-phiy + L*v.y()*phi - fdeltaY*v.z()*phi*std::cos(phi)) - 
     2*fDy1*v.z()*( 2*(fDx4plus2 + fDx3plus1)*fPhiTwist -2*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*phi)*std::sin(phi))/
      (4.*fPhiTwist*v.z()*(4*fDy1*(std::cos(phi) - std::sin(phi)*fTAlph) + fa1md1*std::sin(phi))) ;
    
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
  G4cout << "G4TwistedTrapAlphaSide finished " << G4endl ;
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

G4int G4TwistedTrapAlphaSide::DistanceToSurface(const G4ThreeVector &gp,
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
   G4double uMax = GetValueB(phiR)/2 ;

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

G4int G4TwistedTrapAlphaSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   static const G4double ctol = 0.5 * kCarTolerance;

   G4double phi ;
   G4double yprime ;
   GetPhiUAtX(xx, phi,yprime ) ;

   G4double fYAxisMax =  0.5 * GetValueB(phi) ;
   G4double fYAxisMin =  - fYAxisMax ;

#ifdef G4SPECSDEBUG
   G4cout << "GetAreaCode: phi = " << phi << G4endl ;
   G4cout << "GetAreaCode: yprime = " << yprime << G4endl ;
   G4cout << "Intervall is " << fYAxisMin << " to " << fYAxisMax << G4endl ;
#endif

   G4int areacode = sInside;
   
   if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {

      G4int zaxis = 1;
      
      if (withTol) {

        G4bool isoutside   = false;
        
        // test boundary of yaxis

         if (yprime < fYAxisMin + ctol) {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary; 
            if (yprime <= fYAxisMin - ctol) isoutside = true;

         } else if (yprime > fYAxisMax - ctol) {
            areacode |= (sAxis0 & (sAxisY | sAxisMax)) | sBoundary;
            if (yprime >= fYAxisMax + ctol)  isoutside = true;
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

         if (yprime < fYAxisMin ) {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary;
         } else if (yprime > fYAxisMax) {
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
      G4Exception("G4TwistedTrapAlphaSide::GetAreaCode()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistedTrapAlphaSide::SetCorners()
{

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
    
    G4double x, y, z;
    G4double CosPhi =  std::cos(fPhiTwist/2) ;
    G4double SinPhi =  std::sin(fPhiTwist/2) ;

    // corner of Axis0min and Axis1min

    x = (-( 2*fDy1 + fdeltaY)*SinPhi + CosPhi*( 2*fDx1 - fdeltaX + 2*fDy1*fTAlph))/2. ;
    y = (-( 2*fDy1 + fdeltaY)*CosPhi - SinPhi*( 2*fDx1 - fdeltaX + 2*fDy1*fTAlph))/2. ;
    z = -fDz ;

    SetCorner(sC0Min1Min, x, y, z);
      
    // corner of Axis0max and Axis1min
    x = ( ( 2*fDy1 - fdeltaY)*SinPhi + CosPhi*( 2*fDx2 - fdeltaX - 2*fDy1*fTAlph))/2. ;
    y = ( ( 2*fDy1 - fdeltaY)*CosPhi + SinPhi*(-2*fDx2 + fdeltaX + 2*fDy1*fTAlph))/2. ; 
    z = -fDz ;

    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    x = (-( 2*fDy2 + fdeltaY)*SinPhi + CosPhi*( 2*fDx4 + fdeltaX - 2*fDy2*fTAlph))/2. ;
    y = ( ( 2*fDy2 + fdeltaY)*CosPhi + SinPhi*( 2*fDx4 + fdeltaX - 2*fDy2*fTAlph))/2. ;
    z = fDz ;

    SetCorner(sC0Max1Max, x, y, z);
      
    // corner of Axis0min and Axis1max
    x = ( ( 2*fDy2 - fdeltaY)*SinPhi + CosPhi*( 2*fDx3 + fdeltaX + 2*fDy2*fTAlph))/2. ;
    y = ( (-2*fDy2 + fdeltaY)*CosPhi + SinPhi*( 2*fDx3 + fdeltaX + 2*fDy2*fTAlph))/2. ; 
    z = fDz ;

    SetCorner(sC0Min1Max, x, y, z);

  } else {

    G4Exception("G4TwistedTrapAlphaSide::SetCorners()",
                "NotImplemented", FatalException,
                "Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistedTrapAlphaSide::SetBoundaries()
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
    
  G4Exception("G4TwistedTrapAlphaSide::SetCorners()",
              "NotImplemented", FatalException,
              "Feature NOT implemented !");
  }
  
}



void G4TwistedTrapAlphaSide::GetPhiUAtX( G4ThreeVector p, G4double &phi, G4double &u) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction ( both for a fixed z = pz)
  
  // phi is given by the z coordinate of p

  phi = p.z()/(2*fDz)*fPhiTwist ;

  G4double CosPhi = std::cos(phi) ;
  G4double SinPhi = std::sin(phi) ;
  G4double c0 = 2*(CosPhi - SinPhi*fTAlph) + (fa1md1*SinPhi)/(2*fDy1) ;
  G4double c1 = (-4*fDy1*SinPhi + CosPhi*(fa1md1 - 4*fDy1*fTAlph))/(2*fDy1) ;

  u = ( -(fa1md1*  2*(fDx4plus2  + fDx3plus1 )*fPhiTwist) + 
       2*(fa1md1*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX) - 16*fDy1*fdeltaY)*phi + 
       4*fDy1*(2*(fDx4plus2 + fDx3plus1)*fPhiTwist 
	     - 2*(-2*(fDx4minus2 + fDx3minus1) - 4*fdeltaX)*phi)*fTAlph + 
	8*fPhiTwist*CosPhi*(2*fDx2*p.x() - 2*fDx1*p.x() + 4*fDy1*p.y() - 4*fDy1*p.x()*fTAlph) - 
       8*fPhiTwist*SinPhi*(4*fDy1*p.x() - 2*fDx2*p.y() + 2*fDx1*p.y() + 4*fDy1*p.y()*fTAlph))/
    (8.*fDy1*fPhiTwist*(c0*c0 +  c1*c1))  ;

}


G4ThreeVector G4TwistedTrapAlphaSide::ProjectPoint(const G4ThreeVector &p, 
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


