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
// $Id: G4TwistedTrapBoxSide.cc,v 
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrapBoxSide.cc
//
// Author:
//
//   Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistedTrapBoxSide.hh"
#include "G4JTPolynomialSolver.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistedTrapBoxSide::G4TwistedTrapBoxSide(const G4String     &name,
			   G4double      PhiTwist,    // twist angle
			   G4double      pDz,         // half z lenght
			   G4double      pTheta,      // direction between end planes
			   G4double      pPhi,        // defined by polar and azimutal angles.
			   G4double      pDy1,        // half y length at -pDz
			   G4double      pDx1,        // half x length at -pDz,-pDy
			   G4double      pDy2,        // half y length at +pDz
			   G4double      pDx2,        // half x length at -pDz,+pDy
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
  
  fDy1   = pDy1 ;
  fDy2   = pDy2 ;

  fDz   = pDz ;

  fAlph = pAlph  ;
  fTAlph = std::tan(fAlph) ;

  fTheta = pTheta ;
  fPhi   = pPhi ;

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

G4TwistedTrapBoxSide::~G4TwistedTrapBoxSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistedTrapBoxSide::GetNormal(const G4ThreeVector &tmpxx, 
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

G4int G4TwistedTrapBoxSide::DistanceToSurface(const G4ThreeVector &gp,
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

  G4double phiy    = fPhiTwist * ( v.y()*p.z() - v.z()*p.y() )   ;
  G4double phix    = fPhiTwist * ( v.x()*p.z() - v.z()*p.x() )   ;


  // special case vz = 0

  if ( v.z() == 0. ) {         

    if ( std::fabs(p.z()) <= L ) {     // intersection possible in z
      
      phi = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 

      u = (4*fPhiTwist*(p.y()*v.x() - p.x()*v.y()) + (-4*fdeltaY*v.x()*phi + v.y()*(2*(fDx1 + fDx2)*fPhiTwist + 4*(-fDx1 + fDx2 + fdeltaX)*phi))*
	   std::cos(phi) - (2*(fDx1 + fDx2)*fPhiTwist*v.x() + 4*((-fDx1 + fDx2 + fdeltaX)*v.x() + fdeltaY*v.y())*phi)*std::sin(phi))/
	(4.*fPhiTwist*((v.x() + fTAlph*v.y())*std::cos(phi) + (-(fTAlph*v.x()) + v.y())*std::sin(phi)))   ;

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

    c[7] = 3600*( fPhiTwist*(fDx1 + fDx2)*v.z() + 2*fTAlph*phiy+ 2*phix  ) ;
    c[6] = -7200*(phix*fTAlph  - phiy + 2*fDz*(v.x() + fTAlph*v.y()) - (fdeltaX - fDx1 + fDx2 + fdeltaY*fTAlph ) *v.z()) ;
    c[5] = 300*(48*fDz*(fTAlph*v.x() - v.y()) -  10*phix  +  fPhiTwist*(fDx1 + fDx2)*v.z()   - 10*phiy*fTAlph ) ;
    c[4] = 600*(10*fDz*v.x() + fTAlph*phix + 10*fDz*fTAlph*v.y() - phiy + fdeltaX*v.z() + (fDx2-fDx1)*v.z() + fdeltaY*fTAlph*v.z() ) ;
    c[3] = -1200*fDz*fTAlph*v.x() + 12*phix  + 1200*fDz*v.y() + 12*phiy*fTAlph + 6*(fDx1+fDx2)*fPhiTwist*v.z();
    c[2] = -24*fDz*v.x() - 24*fDz*fTAlph*v.y() + 28*phix*fTAlph + 12*fdeltaX*v.z()  + 12*(fDx2-fDx1)*v.z() + 12*fdeltaY*fTAlph*v.z() - 28*phiy ;
    c[1] = -56*fDz*fTAlph*v.x() + 9*phix + 56*fDz*v.y() + 9*phiy*fTAlph ;
    c[0] = -18*fDz*v.x() - 18*fDz*fTAlph*v.y() ;

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

	u   = - (-4*( -phiy  + L*v.y()*phi) + 4*fdeltaY*v.z()*phi*std::cos(phi) + 
	       v.z()*(2*(fDx1 + fDx2)*fPhiTwist + 4*(-fDx1 + fDx2 + fdeltaX)*phi)*std::sin(phi))/
	  (4.*fPhiTwist*v.z()*(std::cos(phi) - fTAlph*std::sin(phi))) ;

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
	G4Exception("G4TwistedTrapBoxSide::DistanceToSurface()",
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
    u   = - (-4*( -phiy  + L*v.y()*phi) + 4*fdeltaY*v.z()*phi*std::cos(phi) + 
	     v.z()*(2*(fDx1 + fDx2)*fPhiTwist + 4*(-fDx1 + fDx2 + fdeltaX)*phi)*std::sin(phi))/
      (4.*fPhiTwist*v.z()*(std::cos(phi) - fTAlph*std::sin(phi))) ;

    
     
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
    u   = - (-4*( -phiy  + L*v.y()*phi) + 4*fdeltaY*v.z()*phi*std::cos(phi) + 
	     v.z()*(2*(fDx1 + fDx2)*fPhiTwist + 4*(-fDx1 + fDx2 + fdeltaX)*phi)*std::sin(phi))/
      (4.*fPhiTwist*v.z()*(std::cos(phi) - fTAlph*std::sin(phi))) ;

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
  G4cout << "G4TwistedTrapBoxSide finished " << G4endl ;
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

G4int G4TwistedTrapBoxSide::DistanceToSurface(const G4ThreeVector &gp,
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
   G4double uMax = GetValueB(phiR) ;

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

G4int G4TwistedTrapBoxSide::GetAreaCode(const G4ThreeVector &xx, 
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
      G4Exception("G4TwistedTrapBoxSide::GetAreaCode()",
                  "NotImplemented", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistedTrapBoxSide::SetCorners()
{

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis) {
    
    G4double x, y, z;
    G4double CosPhi =  std::cos(fPhiTwist/2) ;
    G4double SinPhi =  std::sin(fPhiTwist/2) ;

    // corner of Axis0min and Axis1min

    x = (-fdeltaX/2. + fDx1 + fDy1*fTAlph)*CosPhi - ((fdeltaY + 2*fDy1)*SinPhi)/2. ;
    y = (-((fdeltaY + 2*fDy1)*CosPhi) + (fdeltaX - 2*(fDx1 + fDy1*fTAlph))*SinPhi)/2. ; 
    z = -fDz ;

    SetCorner(sC0Min1Min, x, y, z);

    // corner of Axis0max and Axis1min

    x = (-fdeltaX/2. + fDx1 - fDy1*fTAlph)*CosPhi - ((fdeltaY - 2*fDy1)*SinPhi)/2. ;
    y = (-fdeltaY/2. + fDy1)*CosPhi +  ((fdeltaX - 2*fDx1 + 2*fDy1*fTAlph)*SinPhi)/2. ;
    z = -fDz ;

    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max

    x = (fdeltaX/2. + fDx2 - fDy2*fTAlph)*CosPhi - ((fdeltaY + 2*fDy2)*SinPhi)/ 2. ;
    y = (fdeltaY/2. + fDy2)*CosPhi + (fdeltaX/2. + fDx2 - fDy2*fTAlph)*SinPhi ;
    z = fDz ;

    SetCorner(sC0Max1Max, x, y, z);
      
    // corner of Axis0min and Axis1max

    x = (fdeltaX/2. + fDx2 + fDy2*fTAlph)*CosPhi - ((fdeltaY - 2*fDy2)*SinPhi)/2. ;
    y = ((fdeltaY - 2*fDy2)*CosPhi)/2. + (fdeltaX/2. + fDx2 + fDy2*fTAlph)*SinPhi ; 
    z = fDz ;

    SetCorner(sC0Min1Max, x, y, z);

  } else {

    G4Exception("G4TwistedTrapBoxSide::SetCorners()",
                "NotImplemented", FatalException,
                "Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistedTrapBoxSide::SetBoundaries()
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
    
  G4Exception("G4TwistedTrapBoxSide::SetCorners()",
              "NotImplemented", FatalException,
              "Feature NOT implemented !");
  }
  
}



void G4TwistedTrapBoxSide::GetPhiUAtX( G4ThreeVector p, G4double &phi, G4double &u) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction ( both for a fixed z = pz)
  
  // phi is given by the z coordinate of p

  phi = p.z()/(2*fDz)*fPhiTwist ;

  u = (-4*fdeltaY*phi + fTAlph*(2*(fDx1 + fDx2)*fPhiTwist + 4*(-fDx1 + fDx2 + fdeltaX)*phi) + 
     4*fPhiTwist*(-(fTAlph*p.x()) + p.y())*std::cos(phi) - 4*fPhiTwist*(p.x() + fTAlph*p.y())*std::sin(phi))/
    (4.*fPhiTwist*(1 + fTAlph*fTAlph )) ;


}


G4ThreeVector G4TwistedTrapBoxSide::ProjectPoint(const G4ThreeVector &p, 
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


