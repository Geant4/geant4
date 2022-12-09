//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4TwistTrapAlphaSide implementation
//
// Author: 18/03/2005 - O.Link (Oliver.Link@cern.ch)
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistTrapAlphaSide.hh"
#include "G4PhysicalConstants.hh"
#include "G4JTPolynomialSolver.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistTrapAlphaSide::
G4TwistTrapAlphaSide(const G4String& name,
                     G4double      PhiTwist,  // twist angle
                     G4double      pDz,       // half z lenght
                     G4double      pTheta,    // direction between end planes
                     G4double      pPhi,      // by polar and azimutal angles
                     G4double      pDy1,      // half y length at -pDz
                     G4double      pDx1,      // half x length at -pDz,-pDy
                     G4double      pDx2,      // half x length at -pDz,+pDy
                     G4double      pDy2,      // half y length at +pDz
                     G4double      pDx3,      // half x length at +pDz,-pDy
                     G4double      pDx4,      // half x length at +pDz,+pDy
                     G4double      pAlph,     // tilt angle at +pDz
                     G4double      AngleSide  // parity
                                             )
  : G4VTwistSurface(name)
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

  fPhiTwist = PhiTwist ;    // dphi
  fAngleSide = AngleSide ;  // 0,90,180,270 deg

  fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi);
    // dx in surface equation
  fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi);
    // dy in surface equation
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;
}


//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistTrapAlphaSide::G4TwistTrapAlphaSide( __void__& a )
  : G4VTwistSurface(a), fTheta(0.), fPhi(0.), fDy1(0.), fDx1(0.), fDx2(0.), 
    fDy2(0.), fDx3(0.), fDx4(0.), fDz(0.), fAlph(0.), fTAlph(0.), fPhiTwist(0.), 
    fAngleSide(0.), fDx4plus2(0.), fDx4minus2(0.), fDx3plus1(0.), fDx3minus1(0.), 
    fDy2plus1(0.), fDy2minus1(0.), fa1md1(0.), fa2md2(0.), fdeltaX(0.),
    fdeltaY(0.)
{
}


//=====================================================================
//* destructor --------------------------------------------------------

G4TwistTrapAlphaSide::~G4TwistTrapAlphaSide()
{
}


//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector
G4TwistTrapAlphaSide::GetNormal(const G4ThreeVector& tmpxx, 
                                      G4bool isGlobal) 
{
   // GetNormal returns a normal vector at a surface (or very close
   // to surface) point at tmpxx.
   // If isGlobal=true, it returns the normal in global coordinate.
   //

   G4ThreeVector xx;
   if (isGlobal)
   {
      xx = ComputeLocalPoint(tmpxx);
      if ((xx - fCurrentNormal.p).mag() < 0.5 * kCarTolerance)
      {
         return ComputeGlobalDirection(fCurrentNormal.normal);
      }
   }
   else
   {
      xx = tmpxx;
      if (xx == fCurrentNormal.p)
      {
         return fCurrentNormal.normal;
      }
   }

   G4double phi ;
   G4double u ;

   GetPhiUAtX(xx,phi,u) ;   // phi,u for point xx close to surface 

   G4ThreeVector normal =  NormAng(phi,u) ;  // the normal vector at phi,u

#ifdef G4TWISTDEBUG
   G4cout  << "normal vector = " << normal << G4endl ;
   G4cout << "phi = " << phi << " , u = " << u << G4endl ;
#endif

   if (isGlobal)
   {
      fCurrentNormal.normal = ComputeGlobalDirection(normal.unit());
   }
   else
   {
      fCurrentNormal.normal = normal.unit();
   }

   return fCurrentNormal.normal;
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int
G4TwistTrapAlphaSide::DistanceToSurface(const G4ThreeVector& gp,
                                        const G4ThreeVector& gv,
                                              G4ThreeVector  gxx[],
                                              G4double       distance[],
                                              G4int          areacode[],
                                              G4bool         isvalid[],
                                              EValidate      validate)
{
  static const G4double pihalf = pi/2 ;
  const G4double ctol = 0.5 * kCarTolerance;

  G4bool IsParallel = false ;
  G4bool IsConverged =  false ;

  G4int nxx = 0 ;  // number of physical solutions

  fCurStatWithV.ResetfDone(validate, &gp, &gv);

  if (fCurStatWithV.IsDone())
  {
    for (G4int i=0; i<fCurStatWithV.GetNXX(); ++i)
    {
      gxx[i] = fCurStatWithV.GetXX(i);
      distance[i] = fCurStatWithV.GetDistance(i);
      areacode[i] = fCurStatWithV.GetAreacode(i);
      isvalid[i]  = fCurStatWithV.IsValid(i);
    }
    return fCurStatWithV.GetNXX();
  }
  else  // initialise
  {
    for (G4int j=0; j<G4VSURFACENXX ; ++j)
    {
      distance[j] = kInfinity;
      areacode[j] = sOutside;
      isvalid[j]  = false;
      gxx[j].set(kInfinity, kInfinity, kInfinity);
    }
  }

  G4ThreeVector p = ComputeLocalPoint(gp);
  G4ThreeVector v = ComputeLocalDirection(gv);
  
#ifdef G4TWISTDEBUG
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

  G4double phixz = fPhiTwist * ( p.x() * v.z() - p.z() * v.x() ) ;
  G4double phiyz = fPhiTwist * ( p.y() * v.z() - p.z() * v.y() ) ;


  // special case vz = 0

  if ( v.z() == 0. )
  {
    if ( std::fabs(p.z()) <= L )   // intersection possible in z
    {
      phi = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 
      u = (fDy1*(4*(-(fdeltaY*phi*v.x()) + fPhiTwist*p.y()*v.x()
                    + fdeltaX*phi*v.y() - fPhiTwist*p.x()*v.y())
                 + ((fDx3plus1 + fDx4plus2)*fPhiTwist
                 +  2*(fDx3minus1 + fDx4minus2)*phi)
                     *(v.y()*std::cos(phi) - v.x()*std::sin(phi))))
         /(fPhiTwist*(4*fDy1* v.x() - (fa1md1 + 4*fDy1*fTAlph)*v.y())
                    *std::cos(phi) + fPhiTwist*(fa1md1*v.x()
                       + 4*fDy1*(fTAlph*v.x() + v.y()))*std::sin(phi));
      xbuftmp.phi = phi ;
      xbuftmp.u = u ;
      xbuftmp.areacode = sOutside ;
      xbuftmp.distance = kInfinity ;
      xbuftmp.isvalid = false ;

      xbuf.push_back(xbuftmp) ;  // store it to xbuf
    }
    else   // no intersection possible
    {
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
  else  // general solution for non-zero vz
  {

    G4double c[8],srd[7],si[7] ;  

    c[7] = 57600*
      fDy1*(fa1md1*phiyz + 
            fDy1*(-4*phixz + 4*fTAlph*phiyz
                 + (fDx3plus1 + fDx4plus2)*fPhiTwist*v.z())) ;
    c[6] = -57600*
      fDy1*(4*fDy1*(phiyz + 2*fDz*v.x() + fTAlph*(phixz - 2*fDz*v.y()))
          - 2*fDy1*(2*fdeltaX + fDx3minus1 + fDx4minus2
                  - 2*fdeltaY*fTAlph)*v.z()
          + fa1md1*(phixz - 2*fDz*v.y() + fdeltaY*v.z()));
    c[5] = 4800*
      fDy1*(fa1md1*(-5*phiyz - 24*fDz*v.x() + 12*fdeltaX*v.z()) + 
            fDy1*(20*phixz - 4*(5*fTAlph*phiyz + 24*fDz*fTAlph*v.x()
                   + 24*fDz*v.y()) + (48*fdeltaY + (fDx3plus1 + fDx4plus2)
                                    *fPhiTwist + 48*fdeltaX*fTAlph)*v.z()));
    c[4] = 4800*
      fDy1*(fa1md1*(phixz - 10*fDz*v.y() + 5*fdeltaY*v.z())
          + 2*fDy1*(2*phiyz + 20*fDz*v.x()
                   + (-10*fdeltaX + fDx3minus1 + fDx4minus2)*v.z()
                   + 2*fTAlph*(phixz - 10*fDz*v.y() + 5*fdeltaY*v.z())));
    c[3] = -96*
      fDy1*(-(fa1md1*(phiyz + 100*fDz*v.x() - 50*fdeltaX*v.z()))
          + fDy1*(4*phixz - 400*fDz*v.y()
                  + (200*fdeltaY - (fDx3plus1 + fDx4plus2)*fPhiTwist)*v.z()
                  - 4*fTAlph*(phiyz + 100*fDz*v.x() - 50*fdeltaX*v.z())));
    c[2] = 32*
      fDy1*(4*fDy1*(7*fTAlph*phixz + 7*phiyz - 6*fDz*v.x() + 6*fDz*fTAlph*v.y())
          + 6*fDy1*(2*fdeltaX+fDx3minus1+fDx4minus2-2*fdeltaY*fTAlph)*v.z()
          + fa1md1*(7*phixz + 6*fDz*v.y() - 3*fdeltaY*v.z()));
    c[1] = -8*
      fDy1*(fa1md1*(-9*phiyz - 56*fDz*v.x() + 28*fdeltaX*v.z())
          + 4*fDy1*(9*phixz - 9*fTAlph*phiyz - 56*fDz*fTAlph*v.x()
                  - 56*fDz*v.y() + 28*(fdeltaY + fdeltaX*fTAlph)*v.z()));
    c[0] = 72*
      fDy1*(fa1md1*(2*fDz*v.y() - fdeltaY*v.z())
          + fDy1*(-8*fDz*v.x() + 8*fDz*fTAlph*v.y()
                + 4*fdeltaX*v.z() - 4*fdeltaY*fTAlph*v.z()));

#ifdef G4TWISTDEBUG
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
    G4int num = trapEq.FindRoots(c,7,srd,si);

    for (G4int i = 0 ; i<num ; i++ )   // loop over all math solutions
    {  
      if ( si[i]==0.0 )  // only real solutions
      { 
#ifdef G4TWISTDEBUG
        G4cout << "Solution " << i << " : " << srd[i] << G4endl ;
#endif
        phi = std::fmod(srd[i] , pihalf)  ;
        u   = (fDy1*(4*(phiyz + 2*fDz*phi*v.y() - fdeltaY*phi*v.z())
                     - ((fDx3plus1 + fDx4plus2)*fPhiTwist
                       + 2*(fDx3minus1 + fDx4minus2)*phi)*v.z()*std::sin(phi)))
             /(fPhiTwist*v.z()*(4*fDy1*std::cos(phi)
                             + (fa1md1 + 4*fDy1*fTAlph)*std::sin(phi)));        
        xbuftmp.phi = phi ;
        xbuftmp.u = u ;
        xbuftmp.areacode = sOutside ;
        xbuftmp.distance = kInfinity ;
        xbuftmp.isvalid = false ;
        
        xbuf.push_back(xbuftmp) ;  // store it to xbuf
      
#ifdef G4TWISTDEBUG
        G4cout << "solution " << i << " = " << phi << " , " << u  << G4endl ;
#endif
      }  // end if real solution
    }  // end loop i    
  }    // end general case

  nxx = (G4int)xbuf.size() ;  // save the number of  solutions

  G4ThreeVector xxonsurface  ;    // point on surface
  G4ThreeVector surfacenormal  ;  // normal vector  
  G4double deltaX;  // distance between intersection point and point on surface
  G4double theta;   // angle between track and surfacenormal
  G4double factor;  // a scaling factor
  G4int maxint=30;  // number of iterations

  for ( std::size_t k = 0 ; k<xbuf.size() ; ++k )
  {
#ifdef G4TWISTDEBUG
    G4cout << "Solution " << k << " : " 
           << "reconstructed phiR = " << xbuf[k].phi
           << ", uR = " << xbuf[k].u << G4endl ; 
#endif
    
    phi = xbuf[k].phi ;  // get the stored values for phi and u
    u = xbuf[k].u ;

    IsConverged = false ;   // no convergence at the beginning
    
    for ( G4int i = 1 ; i<maxint ; ++i )
    {
      xxonsurface = SurfacePoint(phi,u) ;
      surfacenormal = NormAng(phi,u) ;

      tmpdist = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, tmpxx); 
      deltaX = ( tmpxx - xxonsurface ).mag() ; 
      theta = std::fabs(std::acos(v*surfacenormal) - pihalf) ;
      if ( theta < 0.001 )
      { 
        factor = 50 ;
        IsParallel = true ;
      }
      else
      {
        factor = 1 ;
      }

#ifdef G4TWISTDEBUG
      G4cout << "Step i = " << i << ", distance = " << tmpdist
             << ", " << deltaX << G4endl ;
      G4cout << "X = " << tmpxx << G4endl ;
#endif
      
      GetPhiUAtX(tmpxx, phi, u) ;
        // the new point xx is accepted and phi/u replaced
      
#ifdef G4TWISTDEBUG
      G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
      if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
      
    }  // end iterative loop (i)
    
    if ( std::fabs(tmpdist)<ctol ) { tmpdist = 0 ; }

#ifdef G4TWISTDEBUG
    G4cout << "refined solution "  << phi << " , " << u  <<  G4endl ;
    G4cout << "distance = " << tmpdist << G4endl ;
    G4cout << "local X = " << tmpxx << G4endl ;
#endif
    
    tmpisvalid = false ;  // init 

    if ( IsConverged )
    {
      if (validate == kValidateWithTol)
      {
        tmpareacode = GetAreaCode(tmpxx);
        if (!IsOutside(tmpareacode))
        {
          if (tmpdist >= 0) tmpisvalid = true;
        }
      }
      else if (validate == kValidateWithoutTol)
      {
        tmpareacode = GetAreaCode(tmpxx, false);
        if (IsInside(tmpareacode))
        {
          if (tmpdist >= 0) { tmpisvalid = true; }
        }
      }
      else  // kDontValidate
      {
        G4Exception("G4TwistTrapAlphaSide::DistanceToSurface()",
                    "GeomSolids0001", FatalException,
                    "Feature NOT implemented !");
      }
    } 
    else
    {
      tmpdist = kInfinity;     // no convergence after 10 steps 
      tmpisvalid = false ;     // solution is not vaild
    }  

    // store the found values
    // 
    xbuf[k].xx = tmpxx ;
    xbuf[k].distance = tmpdist ;
    xbuf[k].areacode = tmpareacode ;
    xbuf[k].isvalid = tmpisvalid ;

  }  // end loop over physical solutions (variable k)

  std::sort(xbuf.begin() , xbuf.end(), DistanceSort ) ;  // sorting

#ifdef G4TWISTDEBUG
  G4cout << G4endl << "list xbuf after sorting : " << G4endl ;
  G4cout << G4endl << G4endl ;
#endif

  // erase identical intersection (within kCarTolerance) 
  //
  xbuf.erase( std::unique(xbuf.begin(), xbuf.end() , EqualIntersection ),
              xbuf.end() );


  // add guesses
  //
  G4int nxxtmp = (G4int)xbuf.size() ;

  if ( nxxtmp<2 || IsParallel  )  // positive end
  {

#ifdef G4TWISTDEBUG
    G4cout << "add guess at +z/2 .. " << G4endl ;
#endif

    phi = fPhiTwist/2 ;
    u   = 0 ;
    
    xbuftmp.phi = phi ;
    xbuftmp.u = u ;
    xbuftmp.areacode = sOutside ;
    xbuftmp.distance = kInfinity ;
    xbuftmp.isvalid = false ;
    
    xbuf.push_back(xbuftmp) ;  // store it to xbuf

#ifdef G4TWISTDEBUG
    G4cout << "add guess at -z/2 .. " << G4endl ;
#endif

    phi = -fPhiTwist/2 ;
    u   = 0 ;
    
    xbuftmp.phi = phi ;
    xbuftmp.u = u ;
    xbuftmp.areacode = sOutside ;
    xbuftmp.distance = kInfinity ;
    xbuftmp.isvalid = false ;
    
    xbuf.push_back(xbuftmp) ;  // store it to xbuf

    for ( std::size_t k = nxxtmp ; k<xbuf.size() ; ++k )
    {

#ifdef G4TWISTDEBUG
      G4cout << "Solution " << k << " : " 
             << "reconstructed phiR = " << xbuf[k].phi
             << ", uR = " << xbuf[k].u << G4endl ; 
#endif
      
      phi = xbuf[k].phi ;  // get the stored values for phi and u
      u   = xbuf[k].u ;

      IsConverged = false ;   // no convergence at the beginning
      
      for ( G4int i = 1 ; i<maxint ; ++i )
      {
        xxonsurface = SurfacePoint(phi,u) ;
        surfacenormal = NormAng(phi,u) ;
        tmpdist = DistanceToPlaneWithV(p, v, xxonsurface, surfacenormal, tmpxx);
        deltaX = ( tmpxx - xxonsurface ).mag();
        theta = std::fabs(std::acos(v*surfacenormal) - pihalf);
        if ( theta < 0.001 )
        { 
          factor = 50 ;    
        }
        else
        {
          factor = 1 ;
        }
        
#ifdef G4TWISTDEBUG
        G4cout << "Step i = " << i << ", distance = " << tmpdist
               << ", " << deltaX << G4endl
               << "X = " << tmpxx << G4endl ;
#endif

        GetPhiUAtX(tmpxx, phi, u) ;
          // the new point xx is accepted and phi/u replaced
      
#ifdef G4TWISTDEBUG
        G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
        if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
      
      }  // end iterative loop (i)
    
      if ( std::fabs(tmpdist)<ctol )  { tmpdist = 0; }

#ifdef G4TWISTDEBUG
      G4cout << "refined solution "  << phi << " , " << u  <<  G4endl ;
      G4cout << "distance = " << tmpdist << G4endl ;
      G4cout << "local X = " << tmpxx << G4endl ;
#endif

      tmpisvalid = false ;  // init 

      if ( IsConverged )
      {
        if (validate == kValidateWithTol)
        {
          tmpareacode = GetAreaCode(tmpxx);
          if (!IsOutside(tmpareacode))
          {
            if (tmpdist >= 0)  { tmpisvalid = true; }
          }
        }
        else if (validate == kValidateWithoutTol)
        {
          tmpareacode = GetAreaCode(tmpxx, false);
          if (IsInside(tmpareacode))
          {
            if (tmpdist >= 0)  { tmpisvalid = true; }
          }
        }
        else  // kDontValidate
        {
          G4Exception("G4TwistedBoxSide::DistanceToSurface()",
                      "GeomSolids0001", FatalException,
                      "Feature NOT implemented !");
        }
      } 
      else
      {
        tmpdist = kInfinity;     // no convergence after 10 steps 
        tmpisvalid = false ;     // solution is not vaild
      }  

      // store the found values
      //
      xbuf[k].xx = tmpxx ;
      xbuf[k].distance = tmpdist ;
      xbuf[k].areacode = tmpareacode ;
      xbuf[k].isvalid = tmpisvalid ;

    }  // end loop over physical solutions 
  }  // end less than 2 solutions

  // sort again
  std::sort(xbuf.begin() , xbuf.end(), DistanceSort ) ;  // sorting

  // erase identical intersection (within kCarTolerance) 
  xbuf.erase( std::unique(xbuf.begin(), xbuf.end() , EqualIntersection ) ,
              xbuf.end() );

#ifdef G4TWISTDEBUG
  G4cout << G4endl << "list xbuf after sorting : " << G4endl ;
  G4cout << G4endl << G4endl ;
#endif

  nxx = (G4int)xbuf.size() ;   // determine number of solutions again.

  for ( G4int i = 0 ; i<(G4int)xbuf.size() ; ++i )
  {
    distance[i] = xbuf[i].distance;
    gxx[i]      = ComputeGlobalPoint(xbuf[i].xx);
    areacode[i] = xbuf[i].areacode ;
    isvalid[i]  = xbuf[i].isvalid ;
    
    fCurStatWithV.SetCurrentStatus(i, gxx[i], distance[i], areacode[i],
                                     isvalid[i], nxx, validate, &gp, &gv);
#ifdef G4TWISTDEBUG
    G4cout << "element Nr. " << i 
           << ", local Intersection = " << xbuf[i].xx 
           << ", distance = " << xbuf[i].distance 
           << ", u = " << xbuf[i].u 
           << ", phi = " << xbuf[i].phi 
           << ", isvalid = " << xbuf[i].isvalid 
           << G4endl ;
#endif

  }  // end for( i ) loop

#ifdef G4TWISTDEBUG
  G4cout << "G4TwistTrapAlphaSide finished " << G4endl ;
  G4cout << nxx << " possible physical solutions found" << G4endl ;
  for ( G4int k= 0 ; k< nxx ; k++ )
  {
    G4cout << "global intersection Point found: " << gxx[k] << G4endl ;
    G4cout << "distance = " << distance[k] << G4endl ;
    G4cout << "isvalid = " << isvalid[k] << G4endl ;
  }
#endif

  return nxx ;
}


//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int
G4TwistTrapAlphaSide::DistanceToSurface(const G4ThreeVector& gp,
                                              G4ThreeVector  gxx[],
                                              G4double       distance[],
                                              G4int          areacode[])
{  
   const G4double ctol = 0.5 * kCarTolerance;

   fCurStat.ResetfDone(kDontValidate, &gp);

   if (fCurStat.IsDone())
   {
      for (G4int i=0; i<fCurStat.GetNXX(); ++i)
      {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   }
   else  // initialize
   {
      for (G4int i=0; i<G4VSURFACENXX; ++i)
      {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }

   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector xx;  // intersection point
   G4ThreeVector xxonsurface ; // interpolated intersection point 

   // the surfacenormal at that surface point
   //
   G4double phiR = 0  ;
   G4double uR = 0 ;

   G4ThreeVector surfacenormal ; 
   G4double deltaX, uMax ;
   G4double halfphi = 0.5*fPhiTwist ;
   
   for ( G4int i = 1 ; i<20 ; ++i )
   {
     xxonsurface = SurfacePoint(phiR,uR) ;
     surfacenormal = NormAng(phiR,uR) ;
     distance[0] = DistanceToPlane(p,xxonsurface,surfacenormal,xx); // new XX
     deltaX = ( xx - xxonsurface ).mag() ; 

#ifdef G4TWISTDEBUG
     G4cout << "i = " << i << ", distance = " << distance[0]
            << ", " << deltaX << G4endl
            << "X = " << xx << G4endl ;
#endif

     // the new point xx is accepted and phi/psi replaced
     //
     GetPhiUAtX(xx, phiR, uR) ;
     
     if ( deltaX <= ctol ) { break ; }
   }

   // check validity of solution ( valid phi,psi ) 

   uMax = GetBoundaryMax(phiR) ;

   if (  phiR > halfphi ) { phiR =  halfphi ; }
   if ( phiR < -halfphi ) { phiR = -halfphi ; }
   if ( uR > uMax )  { uR = uMax ;  }
   if ( uR < -uMax ) { uR = -uMax ; }

   xxonsurface = SurfacePoint(phiR,uR) ;
   distance[0] = (  p - xx ).mag() ;
   if ( distance[0] <= ctol ) { distance[0] = 0 ; } 

   // end of validity 

#ifdef G4TWISTDEBUG
   G4cout << "refined solution "  << phiR << " , " << uR << " , " << G4endl ;
   G4cout << "distance = " << distance[0] << G4endl ;
   G4cout << "X = " << xx << G4endl ;
#endif

   G4bool isvalid = true;
   gxx[0]      = ComputeGlobalPoint(xx);
   
#ifdef G4TWISTDEBUG
   G4cout << "intersection Point found: " << gxx[0] << G4endl ;
   G4cout << "distance = " << distance[0] << G4endl ;
#endif

   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                            isvalid, 1, kDontValidate, &gp);
   return 1;
}


//=====================================================================
//* GetAreaCode -------------------------------------------------------

G4int
G4TwistTrapAlphaSide::GetAreaCode(const G4ThreeVector& xx, G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   const G4double ctol = 0.5 * kCarTolerance;

   G4double phi ;
   G4double yprime ;
   GetPhiUAtX(xx, phi,yprime ) ;

   G4double fYAxisMax =  GetBoundaryMax(phi) ;
   G4double fYAxisMin =  GetBoundaryMin(phi) ;

#ifdef G4TWISTDEBUG
   G4cout << "GetAreaCode: phi = " << phi << G4endl ;
   G4cout << "GetAreaCode: yprime = " << yprime << G4endl ;
   G4cout << "Intervall is " << fYAxisMin << " to " << fYAxisMax << G4endl ;
#endif

   G4int areacode = sInside;
   
   if (fAxis[0] == kYAxis && fAxis[1] == kZAxis)
   {
      G4int zaxis = 1;

      if (withTol)
      {
        G4bool isoutside   = false;
        
        // test boundary of yaxis

         if (yprime < fYAxisMin + ctol)
         {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary; 
            if (yprime <= fYAxisMin - ctol)  { isoutside = true; }

         }
         else if (yprime > fYAxisMax - ctol)
         {
            areacode |= (sAxis0 & (sAxisY | sAxisMax)) | sBoundary;
            if (yprime >= fYAxisMax + ctol)  { isoutside = true; }
         }

         // test boundary of z-axis

         if (xx.z() < fAxisMin[zaxis] + ctol)
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin)); 

            if   (areacode & sBoundary)   // xx is on the corner
              { areacode |= sCorner; }

            else
              { areacode |= sBoundary; }
            if (xx.z() <= fAxisMin[zaxis] - ctol)  { isoutside = true; }
         }
         else if (xx.z() > fAxisMax[zaxis] - ctol)
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax));

            if   (areacode & sBoundary)   // xx is on the corner
               { areacode |= sCorner; }
            else
               { areacode |= sBoundary; }
            if (xx.z() >= fAxisMax[zaxis] + ctol)  { isoutside = true; }
         }

         // if isoutside = true, clear inside bit.             
         // if not on boundary, add axis information.             
         
         if (isoutside)
         {
            G4int tmpareacode = areacode & (~sInside);
            areacode = tmpareacode;
         }
         else if ((areacode & sBoundary) != sBoundary)
         {
            areacode |= (sAxis0 & sAxisY) | (sAxis1 & sAxisZ);
         }           
         
      }
      else
      {
         // boundary of y-axis

         if (yprime < fYAxisMin )
         {
            areacode |= (sAxis0 & (sAxisY | sAxisMin)) | sBoundary;
         }
         else if (yprime > fYAxisMax)
         {
            areacode |= (sAxis0 & (sAxisY | sAxisMax)) | sBoundary;
         }
         
         // boundary of z-axis

         if (xx.z() < fAxisMin[zaxis])
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin));
            if   (areacode & sBoundary)   // xx is on the corner
              { areacode |= sCorner; }
            else
              { areacode |= sBoundary; }
         }
         else if (xx.z() > fAxisMax[zaxis])
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax)) ;
            if   (areacode & sBoundary)   // xx is on the corner
              { areacode |= sCorner; }
            else
              { areacode |= sBoundary; }
         }

         if ((areacode & sBoundary) != sBoundary)
         {
            areacode |= (sAxis0 & sAxisY) | (sAxis1 & sAxisZ);
         }           
      }
      return areacode;
   }
   else
   {
      G4Exception("G4TwistTrapAlphaSide::GetAreaCode()",
                  "GeomSolids0001", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistTrapAlphaSide::SetCorners()
{

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis)
  {
    
    G4double x, y, z;

    // corner of Axis0min and Axis1min
    //
    x = -fdeltaX/2. + (fDx1 - fDy1*fTAlph)*std::cos(fPhiTwist/2.)
      - fDy1*std::sin(fPhiTwist/2.);
    y = -fdeltaY/2. - fDy1*std::cos(fPhiTwist/2.)
      + (-fDx1 + fDy1*fTAlph)*std::sin(fPhiTwist/2.);
    z = -fDz ;

    //    G4cout << "SetCorners: " << x << ", " << y << ", " << z << G4endl ;

    SetCorner(sC0Min1Min, x, y, z);
      
    // corner of Axis0max and Axis1min
    //
    x = -fdeltaX/2. + (fDx2 + fDy1*fTAlph)*std::cos(fPhiTwist/2.)
      + fDy1*std::sin(fPhiTwist/2.);
    y = -fdeltaY/2. + fDy1*std::cos(fPhiTwist/2.)
      - (fDx2 + fDy1*fTAlph)*std::sin(fPhiTwist/2.);
    z = -fDz ;

    //    G4cout << "SetCorners: " << x << ", " << y << ", " << z << G4endl ;

    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    //
    x = fdeltaX/2. + (fDx4 + fDy2*fTAlph)*std::cos(fPhiTwist/2.)
      - fDy2*std::sin(fPhiTwist/2.);
    y = fdeltaY/2. + fDy2*std::cos(fPhiTwist/2.)
      + (fDx4 + fDy2*fTAlph)*std::sin(fPhiTwist/2.);
    z = fDz ;
    
    //    G4cout << "SetCorners: " << x << ", " << y << ", " << z << G4endl ;
    
    SetCorner(sC0Max1Max, x, y, z);

    // corner of Axis0min and Axis1max
    x = fdeltaX/2. + (fDx3 - fDy2*fTAlph)*std::cos(fPhiTwist/2.)
      + fDy2*std::sin(fPhiTwist/2.) ;
    y = fdeltaY/2. - fDy2*std::cos(fPhiTwist/2.)
      + (fDx3 - fDy2*fTAlph)*std::sin(fPhiTwist/2.) ;
    z = fDz ;

    //    G4cout << "SetCorners: " << x << ", " << y << ", " << z << G4endl ;

    SetCorner(sC0Min1Max, x, y, z);

  }
  else
  {
    G4Exception("G4TwistTrapAlphaSide::SetCorners()",
                "GeomSolids0001", FatalException,
                "Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistTrapAlphaSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   

  G4ThreeVector direction;
   
  if (fAxis[0] == kYAxis && fAxis[1] == kZAxis)
  {     
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
    
  }
  else
  {
    G4Exception("G4TwistTrapAlphaSide::SetCorners()",
                "GeomSolids0001", FatalException,
                "Feature NOT implemented !");
  }
}

//=====================================================================
//* GetPhiUAtX --------------------------------------------------------

void
G4TwistTrapAlphaSide::GetPhiUAtX( G4ThreeVector p, G4double& phi, G4double& u ) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction
  // ( both for a fixed z = pz)
  
  // phi is given by the z coordinate of p

  phi =  p.z()/(2*fDz)*fPhiTwist ;
  u = (fPhiTwist*(2*fDx1*fDx1 - 2*fDx2*fDx2 - fa1md1*(fDx3 + fDx4)
             - 4*(fDx3plus1 + fDx4plus2)*fDy1*fTAlph)
             - 2*(2*fDx1*fDx1 - 2*fDx2*fDx2 + fa1md1*(fDx3 + fDx4)
                + 4*(fDx3minus1 + fDx4minus2)*fDy1*fTAlph)*phi
             - 4*(fa1md1*(fdeltaX*phi - fPhiTwist*p.x())
             + 4*fDy1*(fdeltaY*phi + fdeltaX*fTAlph*phi
                     - fPhiTwist*(fTAlph*p.x() + p.y())))*std::cos(phi)
             - 4*(fa1md1*fdeltaY*phi - 4*fdeltaX*fDy1*phi
                + 4*fdeltaY*fDy1*fTAlph*phi + 4*fDy1*fPhiTwist*p.x()
                - fPhiTwist*(fa1md1 + 4*fDy1*fTAlph)*p.y())*std::sin(phi))
       /(fDy1* fPhiTwist*((std::fabs(((fa1md1 + 4*fDy1*fTAlph)*std::cos(phi))
                          /fDy1 - 4*std::sin(phi)))
                        *(std::fabs(((fa1md1 + 4*fDy1*fTAlph)*std::cos(phi))
                          /fDy1 - 4*std::sin(phi)))
             + (std::fabs(4*std::cos(phi)
                + ((fa1md1 + 4*fDy1*fTAlph)*std::sin(phi))/fDy1))
             * (std::fabs(4*std::cos(phi)
                + ((fa1md1 + 4*fDy1*fTAlph)*std::sin(phi))/fDy1)))) ;
}

//=====================================================================
//* ProjectPoint ------------------------------------------------------

G4ThreeVector
G4TwistTrapAlphaSide::ProjectPoint(const G4ThreeVector& p, G4bool isglobal) 
{
  // Get Rho at p.z() on Hyperbolic Surface.

  G4ThreeVector tmpp;
  if (isglobal)
  {
     tmpp = fRot.inverse()*p - fTrans;
  }
  else
  {
     tmpp = p;
  }

  G4double phi ;
  G4double u ;

  GetPhiUAtX( tmpp, phi, u ) ;
    // calculate (phi, u) for a point p close the surface
  
  G4ThreeVector xx = SurfacePoint(phi,u) ;
    // transform back to cartesian coordinates

  if (isglobal)
  {
     return (fRot * xx + fTrans);
  }
  else
  {
     return xx;
  }
}

//=====================================================================
//* GetFacets ---------------------------------------------------------

void
G4TwistTrapAlphaSide::GetFacets( G4int k, G4int n, G4double xyz[][3],
                                 G4int faces[][4], G4int iside ) 
{

  G4double phi ;
  G4double b ;   

  G4double z, u ;     // the two parameters for the surface equation
  G4ThreeVector p ;  // a point on the surface, given by (z,u)

  G4int nnode ;
  G4int nface ;

  // calculate the (n-1)*(k-1) vertices

  for ( G4int i = 0 ; i<n ; ++i )
  {
    z = -fDz+i*(2.*fDz)/(n-1) ;
    phi = z*fPhiTwist/(2*fDz) ;
    b = GetValueB(phi) ;

    for ( G4int j = 0 ; j<k ; ++j )
    {
      nnode = GetNode(i,j,k,n,iside) ;
      u = -b/2 +j*b/(k-1) ;
      p = SurfacePoint(phi,u,true) ;  // surface point in global coordinates

      xyz[nnode][0] = p.x() ;
      xyz[nnode][1] = p.y() ;
      xyz[nnode][2] = p.z() ;

      if ( i<n-1 && j<k-1 )  // conterclock wise filling 
      {
        nface = GetFace(i,j,k,n,iside) ;
        faces[nface][0] = GetEdgeVisibility(i,j,k,n,0,-1)
                        * (GetNode(i  ,j  ,k,n,iside)+1) ;  // f77 numbering
        faces[nface][1] = GetEdgeVisibility(i,j,k,n,1,-1)
                        * (GetNode(i  ,j+1,k,n,iside)+1) ;
        faces[nface][2] = GetEdgeVisibility(i,j,k,n,2,-1)
                        * (GetNode(i+1,j+1,k,n,iside)+1) ;
        faces[nface][3] = GetEdgeVisibility(i,j,k,n,3,-1)
                        * (GetNode(i+1,j  ,k,n,iside)+1) ;
      }
    }
  }
}
