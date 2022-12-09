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
// G4TwistTrapParallelSide implementation
//
// Author: Oliver Link (Oliver.Link@cern.ch)
// --------------------------------------------------------------------

#include <cmath>

#include "G4TwistTrapParallelSide.hh"
#include "G4PhysicalConstants.hh"
#include "G4JTPolynomialSolver.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistTrapParallelSide::G4TwistTrapParallelSide(const G4String& name,
                           G4double PhiTwist,  // twist angle
                           G4double pDz,       // half z lenght
                           G4double pTheta,    // direction between end planes
                           G4double pPhi,      // by polar and azimutal angles
                           G4double pDy1,      // half y length at -pDz
                           G4double pDx1,      // half x length at -pDz,-pDy
                           G4double pDx2,      // half x length at -pDz,+pDy
                           G4double pDy2,      // half y length at +pDz
                           G4double pDx3,      // half x length at +pDz,-pDy
                           G4double pDx4,      // half x length at +pDz,+pDy
                           G4double pAlph,     // tilt angle at +pDz
                           G4double AngleSide  // parity
                                               )
  : G4VTwistSurface(name)
{  
  
  fAxis[0]    = kXAxis; // in local coordinate system
  fAxis[1]    = kZAxis;
  fAxisMin[0] = -kInfinity ;  // X Axis boundary
  fAxisMax[0] = kInfinity ;   //   depends on z !!
  fAxisMin[1] = -pDz ;        // Z Axis boundary
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
  //
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

  fdeltaX = 2*fDz*std::tan(fTheta)*std::cos(fPhi); // dx in surface equation
  fdeltaY = 2*fDz*std::tan(fTheta)*std::sin(fPhi); // dy in surface equation
  
  fRot.rotateZ( AngleSide ) ; 
  
  fTrans.set(0, 0, 0);  // No Translation
  fIsValidNorm = false;
  
  SetCorners() ;
  SetBoundaries() ;
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistTrapParallelSide::G4TwistTrapParallelSide( __void__& a )
  : G4VTwistSurface(a), fTheta(0.), fPhi(0.), fDy1(0.), fDx1(0.), fDx2(0.), 
    fDy2(0.), fDx3(0.), fDx4(0.), fDz(0.), fAlph(0.), fTAlph(0.), fPhiTwist(0.),
    fAngleSide(0.), fdeltaX(0.), fdeltaY(0.), fDx4plus2(0.), fDx4minus2(0.), 
    fDx3plus1(0.), fDx3minus1(0.), fDy2plus1(0.), fDy2minus1(0.), fa1md1(0.), 
    fa2md2(0.)
{
}

//=====================================================================
//* destructor --------------------------------------------------------

G4TwistTrapParallelSide::~G4TwistTrapParallelSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistTrapParallelSide::GetNormal(const G4ThreeVector& tmpxx, 
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

   //    normal = normal/normal.mag() ;

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

G4int G4TwistTrapParallelSide::DistanceToSurface(const G4ThreeVector& gp,
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
  else   // initialize
  {
    for (G4int i=0; i<G4VSURFACENXX ; ++i)
    {
      distance[i] = kInfinity;
      areacode[i] = sOutside;
      isvalid[i]  = false;
      gxx[i].set(kInfinity, kInfinity, kInfinity);
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
    if ( std::fabs(p.z()) <= L )       // intersection possible in z
    {
      phi = p.z() * fPhiTwist / L ;  // phi is determined by the z-position 

      u = (2*(fdeltaY*phi*v.x() - fPhiTwist*p.y()*v.x() - fdeltaX*phi*v.y()
              + fPhiTwist*p.x()*v.y()) + (fDy2plus1*fPhiTwist
              + 2*fDy2minus1*phi)*(v.x()*std::cos(phi) + v.y()*std::sin(phi)))
        / (2.* fPhiTwist*(v.y()*std::cos(phi) - v.x()*std::sin(phi)));

      xbuftmp.phi = phi ;
      xbuftmp.u = u ;
      xbuftmp.areacode = sOutside ;
      xbuftmp.distance = kInfinity ;
      xbuftmp.isvalid = false ;

      xbuf.push_back(xbuftmp) ; // store it to xbuf
    }
    else                              // no intersection possible
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
  else   // general solution for non-zero vz
  {
    G4double c[9],srd[8],si[8] ;  

    c[8] = -3600*(-2*phiyz + fDy2plus1*fPhiTwist*v.z()) ;
    c[7] = -7200*(phixz - 2*fDz*v.y() + (fdeltaY + fDy2minus1)*v.z()) ;
    c[6] = 120*(-52*phiyz - 120*fDz*v.x() + 60*fdeltaX*v.z()
         + 11*fDy2plus1*fPhiTwist*v.z()) ;
    c[5] = 240*(16*phixz - 52*fDz*v.y() + 26*fdeltaY*v.z()
         + 11*fDy2minus1*v.z()) ;
    c[4] = 12*(127*phiyz + 640*fDz*v.x() - 320*fdeltaX*v.z()
         + 4*fDy2plus1*fPhiTwist*v.z()) ;
    c[3] = -404*phixz + 3048*fDz*v.y() - 1524*fdeltaY*v.z()
         + 96*fDy2minus1*v.z() ;
    c[2] = -72*phiyz + 404*(-2*fDz*v.x() + fdeltaX*v.z()) ;
    c[1] = 12*(phixz - 12*fDz*v.y() + 6*fdeltaY*v.z()) ;
    c[0] = 24*fDz*v.x() - 12*fdeltaX*v.z() ;


#ifdef G4TWISTDEBUG
    G4cout << "coef = " << c[0] << " " 
           <<  c[1] << " "  
           <<  c[2] << " "  
           <<  c[3] << " "  
           <<  c[4] << " "  
           <<  c[5] << " "  
           <<  c[6] << " "  
           <<  c[7] << " "  
           <<  c[8] << G4endl ;
#endif    

    G4JTPolynomialSolver trapEq ;
    G4int num = trapEq.FindRoots(c,8,srd,si);

    for (G4int i = 0 ; i<num ; ++i )   // loop over all mathematical solutions
    {
      if ( si[i]==0.0 )   // only real solutions
      {
#ifdef G4TWISTDEBUG
        G4cout << "Solution " << i << " : " << srd[i] << G4endl ;
#endif
        phi = std::fmod(srd[i] , pihalf)  ;
        u = (1/std::cos(phi)*(2*phixz + 4*fDz*phi*v.x()
          - 2*fdeltaX*phi*v.z() + (fDy2plus1*fPhiTwist
          + 2*fDy2minus1*phi)*v.z()* std::sin(phi)))/(2.*fPhiTwist*v.z()) ;

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
  }  // end general case

  nxx = (G4int)xbuf.size() ;  // save the number of  solutions

  G4ThreeVector xxonsurface  ;       // point on surface
  G4ThreeVector surfacenormal  ;     // normal vector  
  G4double deltaX  ; // distance between intersection point and point on surface
  G4double theta  ;                  // angle between track and surfacenormal
  G4double factor ;                  // a scaling factor
  G4int maxint = 30 ;                // number of iterations


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
      G4cout << "Step i = " << i << ", distance = "
             << tmpdist << ", " << deltaX << G4endl ;
      G4cout << "X = " << tmpxx << G4endl ;
#endif
      
      GetPhiUAtX(tmpxx, phi, u); // new point xx is accepted and phi/u replaced
      
#ifdef G4TWISTDEBUG
      G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
      if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
    }  // end iterative loop (i)
    
    if ( std::fabs(tmpdist)<ctol ) tmpdist = 0 ; 

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
          if (tmpdist >= 0) tmpisvalid = true;
        }
      }
      else  // kDontValidate
      {
        G4Exception("G4TwistTrapParallelSide::DistanceToSurface()",
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
  xbuf.erase(std::unique(xbuf.begin(),xbuf.end(),EqualIntersection),xbuf.end());


  // add guesses

  G4int nxxtmp = (G4int)xbuf.size() ;

  if ( nxxtmp<2 || IsParallel  )
  {
    // positive end
#ifdef G4TWISTDEBUG
    G4cout << "add guess at +z/2 .. " << G4endl ;
#endif

    phi = fPhiTwist/2 ;
    u   =  0 ;
    
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
        deltaX = ( tmpxx - xxonsurface ).mag() ; 
        theta = std::fabs(std::acos(v*surfacenormal) - pihalf) ;
        if ( theta < 0.001 )
        { 
          factor = 50 ;    
        }
        else
        {
          factor = 1 ;
        }
        
#ifdef G4TWISTDEBUG
        G4cout << "Step i = " << i << ", distance = "
               << tmpdist << ", " << deltaX << G4endl ;
        G4cout << "X = " << tmpxx << G4endl ;
#endif

        GetPhiUAtX(tmpxx, phi, u) ; // new point xx accepted and phi/u replaced

#ifdef G4TWISTDEBUG
        G4cout << "approximated phi = " << phi << ", u = " << u << G4endl ; 
#endif
      
        if ( deltaX <= factor*ctol ) { IsConverged = true ; break ; }
      }  // end iterative loop (i)

      if ( std::fabs(tmpdist)<ctol ) tmpdist = 0 ; 

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
            if (tmpdist >= 0) tmpisvalid = true;
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
      xbuf[k].xx = tmpxx ;
      xbuf[k].distance = tmpdist ;
      xbuf[k].areacode = tmpareacode ;
      xbuf[k].isvalid = tmpisvalid ;

    }  // end loop over physical solutions 
  }  // end less than 2 solutions

  // sort again
  std::sort(xbuf.begin() , xbuf.end(), DistanceSort ) ;  // sorting

  // erase identical intersection (within kCarTolerance) 
  xbuf.erase(std::unique(xbuf.begin(),xbuf.end(),EqualIntersection),xbuf.end());

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
  G4cout << "G4TwistTrapParallelSide finished " << G4endl ;
  G4cout << nxx << " possible physical solutions found" << G4endl ;
  for ( G4int k= 0 ; k< nxx ; ++k )
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

G4int G4TwistTrapParallelSide::DistanceToSurface(const G4ThreeVector& gp,
                                                 G4ThreeVector gxx[],
                                                 G4double      distance[],
                                                 G4int         areacode[])
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
   G4double phiR = 0  ; // 
   G4double uR = 0 ;

   G4ThreeVector surfacenormal ; 
   G4double deltaX ;
   
   G4int maxint = 20 ;

   for ( G4int i = 1 ; i<maxint ; ++i )
   {
     xxonsurface = SurfacePoint(phiR,uR) ;
     surfacenormal = NormAng(phiR,uR) ;
     distance[0] = DistanceToPlane(p, xxonsurface, surfacenormal, xx); // new XX
     deltaX = ( xx - xxonsurface ).mag() ; 

#ifdef G4TWISTDEBUG
     G4cout << "i = " << i << ", distance = "
            << distance[0] << ", " << deltaX << G4endl ;
     G4cout << "X = " << xx << G4endl ;
#endif

     // the new point xx is accepted and phi/psi replaced
     GetPhiUAtX(xx, phiR, uR) ;
     
     if ( deltaX <= ctol ) { break ; }
   }

   // check validity of solution ( valid phi,psi ) 

   G4double halfphi = 0.5*fPhiTwist ;
   G4double uMax = GetBoundaryMax(phiR) ;
   G4double uMin = GetBoundaryMin(phiR) ;

   if (  phiR > halfphi ) phiR =  halfphi ;
   if ( phiR < -halfphi ) phiR = -halfphi ;
   if ( uR > uMax ) uR = uMax ;
   if ( uR < uMin ) uR = uMin ;

   xxonsurface = SurfacePoint(phiR,uR) ;
   distance[0] = (  p - xx ).mag() ;
   if ( distance[0] <= ctol ) { distance[0] = 0 ; } 

   // end of validity 

#ifdef G4TWISTDEBUG
   G4cout << "refined solution "  << phiR << " , " << uR << " , " <<  G4endl ;
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

G4int G4TwistTrapParallelSide::GetAreaCode(const G4ThreeVector& xx, 
                                           G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   const G4double ctol = 0.5 * kCarTolerance;

   G4double phi ;
   G4double yprime ;
   GetPhiUAtX(xx, phi,yprime ) ;

   G4double fXAxisMax = GetBoundaryMax(phi) ;
   G4double fXAxisMin = GetBoundaryMin(phi) ;

#ifdef G4TWISTDEBUG
   G4cout << "GetAreaCode: phi = " << phi << G4endl ;
   G4cout << "GetAreaCode: yprime = " << yprime << G4endl ;
   G4cout << "Intervall is " << fXAxisMin << " to " << fXAxisMax << G4endl ;
#endif

   G4int areacode = sInside;
   
   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis)
   {
      G4int zaxis = 1;
      
      if (withTol)
      {
        G4bool isoutside = false;
        
        // test boundary of xaxis

         if (yprime < fXAxisMin + ctol)
         {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary; 
            if (yprime <= fXAxisMin - ctol) isoutside = true;

         }
         else if (yprime > fXAxisMax - ctol)
         {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
            if (yprime >= fXAxisMax + ctol)  isoutside = true;
         }

         // test boundary of z-axis

         if (xx.z() < fAxisMin[zaxis] + ctol)
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin)); 

            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on corner.
            else                        areacode |= sBoundary;
            if (xx.z() <= fAxisMin[zaxis] - ctol) isoutside = true;

         }
         else if (xx.z() > fAxisMax[zaxis] - ctol)
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax));

            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on corner.
            else                        areacode |= sBoundary; 
            if (xx.z() >= fAxisMax[zaxis] + ctol) isoutside = true;
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
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
         
      }
      else
      {
         // boundary of y-axis

         if (yprime < fXAxisMin )
         {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary;
         }
         else if (yprime > fXAxisMax)
         {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
         }
         
         // boundary of z-axis

         if (xx.z() < fAxisMin[zaxis])
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on corner.
            else                        areacode |= sBoundary; 
           
         }
         else if (xx.z() > fAxisMax[zaxis])
         {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax)) ;
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on corner.
            else                        areacode |= sBoundary; 
         }

         if ((areacode & sBoundary) != sBoundary)
         {
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
      }
      return areacode;
   }
   else
   {
      G4Exception("G4TwistTrapParallelSide::GetAreaCode()",
                  "GeomSolids0001", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistTrapParallelSide::SetCorners()
{

  // Set Corner points in local coodinate.   

  if (fAxis[0] == kXAxis && fAxis[1] == kZAxis)
  {
    G4double x, y, z;

    // corner of Axis0min and Axis1min

    x = -fdeltaX/2. + (-fDx2 + fDy1*fTAlph)*std::cos(fPhiTwist/2.)
      + fDy1*std::sin(fPhiTwist/2.) ;
    y = -fdeltaY/2. + fDy1*std::cos(fPhiTwist/2.)
      + (fDx2 - fDy1*fTAlph)*std::sin(fPhiTwist/2.) ;
    z = -fDz ;

    SetCorner(sC0Min1Min, x, y, z);
      
    // corner of Axis0max and Axis1min

    x = -fdeltaX/2. + (fDx2 + fDy1*fTAlph)*std::cos(fPhiTwist/2.)
      + fDy1*std::sin(fPhiTwist/2.) ;
    y = -fdeltaY/2. + fDy1*std::cos(fPhiTwist/2.)
      - (fDx2 + fDy1*fTAlph)*std::sin(fPhiTwist/2.) ;
    z = -fDz;

    SetCorner(sC0Max1Min, x, y, z);
      
    // corner of Axis0max and Axis1max
    x = fdeltaX/2. + (fDx4 + fDy2*fTAlph)*std::cos(fPhiTwist/2.)
      - fDy2*std::sin(fPhiTwist/2.) ;
    y = fdeltaY/2. + fDy2*std::cos(fPhiTwist/2.)
      + (fDx4 + fDy2*fTAlph)*std::sin(fPhiTwist/2.) ;
    z = fDz ;

    SetCorner(sC0Max1Max, x, y, z);
      
    // corner of Axis0min and Axis1max
    x = fdeltaX/2. + (-fDx4 + fDy2*fTAlph)*std::cos(fPhiTwist/2.)
      - fDy2*std::sin(fPhiTwist/2.) ;
    y = fdeltaY/2. + fDy2*std::cos(fPhiTwist/2.)
      + (-fDx4 + fDy2*fTAlph)*std::sin(fPhiTwist/2.) ;
    z = fDz ;

    SetCorner(sC0Min1Max, x, y, z);
  }
  else
  {
    G4Exception("G4TwistTrapParallelSide::SetCorners()",
                "GeomSolids0001", FatalException,
                "Method NOT implemented !");
  }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistTrapParallelSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   

  G4ThreeVector direction;
   
  if (fAxis[0] == kXAxis && fAxis[1] == kZAxis)
  {
    // sAxis0 & sAxisMin
    direction = GetCorner(sC0Min1Max) - GetCorner(sC0Min1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisX | sAxisMin), direction, 
                GetCorner(sC0Min1Min), sAxisZ) ;
      
      // sAxis0 & sAxisMax
    direction = GetCorner(sC0Max1Max) - GetCorner(sC0Max1Min);
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisX | sAxisMax), direction, 
                GetCorner(sC0Max1Min), sAxisZ);
    
    // sAxis1 & sAxisMin
    direction = GetCorner(sC0Max1Min) - GetCorner(sC0Min1Min);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
                GetCorner(sC0Min1Min), sAxisX);
    
    // sAxis1 & sAxisMax
    direction = GetCorner(sC0Max1Max) - GetCorner(sC0Min1Max);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
                GetCorner(sC0Min1Max), sAxisX);
  }
  else
  {
    G4Exception("G4TwistTrapParallelSide::SetCorners()",
                "GeomSolids0001", FatalException,
                "Feature NOT implemented !");
  }
}

//=====================================================================
//* GetPhiUAtX() ------------------------------------------------------

void
G4TwistTrapParallelSide::GetPhiUAtX( G4ThreeVector p,
                                     G4double& phi, G4double& u ) 
{
  // find closest point XX on surface for a given point p
  // X0 is a point on the surface,  d is the direction
  // ( both for a fixed z = pz)
  
  // phi is given by the z coordinate of p

  phi = p.z()/(2*fDz)*fPhiTwist ;

  u = ((-(fdeltaX*phi) + fPhiTwist*p.x())* std::cos(phi)
    + (-(fdeltaY*phi) + fPhiTwist*p.y())*std::sin(phi))/fPhiTwist ;
}

//=====================================================================
//* ProjectPoint() ----------------------------------------------------

G4ThreeVector G4TwistTrapParallelSide::ProjectPoint(const G4ThreeVector& p, 
                                                    G4bool isglobal) 
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

  GetPhiUAtX( tmpp, phi, u ) ;  // calculate (phi, u) for p close to surface
  
  G4ThreeVector xx = SurfacePoint(phi,u) ; // transform back to Cartesian coords

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
//* GetFacets() -------------------------------------------------------

void G4TwistTrapParallelSide::GetFacets( G4int k, G4int n, G4double xyz[][3],
                                         G4int faces[][4], G4int iside ) 
{
  G4double phi ;
  G4double z, u ;     // the two parameters for the surface equation
  G4ThreeVector p ;  // a point on the surface, given by (z,u)

  G4int nnode ;
  G4int nface ;

  G4double umin, umax ;

  // calculate the (n-1)*(k-1) vertices

  for ( G4int i = 0 ; i<n ; ++i )
  {
    z = -fDz+i*(2.*fDz)/(n-1) ;
    phi = z*fPhiTwist/(2*fDz) ;
    umin = GetBoundaryMin(phi) ;
    umax = GetBoundaryMax(phi) ;

    for ( G4int j = 0 ; j<k ; ++j )
    {
      nnode = GetNode(i,j,k,n,iside) ;
      u = umax - j*(umax-umin)/(k-1) ;
      p = SurfacePoint(phi,u,true) ;  // surface point in global coords

      xyz[nnode][0] = p.x() ;
      xyz[nnode][1] = p.y() ;
      xyz[nnode][2] = p.z() ;

      if ( i<n-1 && j<k-1 )    // conterclock wise filling
      {
        nface = GetFace(i,j,k,n,iside) ;
        faces[nface][0] = GetEdgeVisibility(i,j,k,n,0,-1)
                        * (GetNode(i  ,j  ,k,n,iside)+1) ; // fortran numbering
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
