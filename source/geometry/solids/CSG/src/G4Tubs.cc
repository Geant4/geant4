// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Tubs.cc,v 1.14 2000-04-11 16:04:30 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Tubs
//
// History:
//
// 1994-95  P.Kent,     implementation
//
// 18.06.98 V.Grichine, n-normalisation in DistanceToOut(p.v) 
// 09.10.98 V.Grichine, modifications in Distance ToOut(p,v,...)
// 23.03.99 V.Grichine, bug fixed in DistanceToIn(p,v) 
// 25.05.99 V.Grichine, bugs fixed in DistanceToIn(p,v) 
// 28.05.99 V.Grichine, bugs fixed in  Distance ToOut(p,v,...)
// 13.10.99 V.Grichine, bugs fixed in DistanceToIn(p,v) 
// 19.11.99 V. Grichine, side = kNull in Distance ToOut(p,v,...)
// 31.03.00 V.Grichine, bug fixed in Inside(p)


#include "G4Tubs.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBStube.hh"
#include "G4NURBScylinder.hh"
#include "G4NURBStubesector.hh"

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4Tubs::G4Tubs( const G4String &pName,
	              G4double pRMin, G4double pRMax,
	              G4double pDz,
	              G4double pSPhi, G4double pDPhi)
: G4CSGSolid(pName)
{

  if (pDz>0) // Check z-len
  {
    fDz = pDz ;
  }
  else
  {
    G4Exception("Error in G4Tubs::G4Tubs - invalid z half-length");
  }
  if ( pRMin < pRMax && pRMin >= 0 ) // Check radii
  {
    fRMin = pRMin ; 
    fRMax = pRMax ;
  }
  else
  {
    G4Exception("Error in G4Tubs::G4Tubs - invalid radii");
  }
  if ( pDPhi >= 2.0*M_PI ) // Check angles
  {
    fDPhi=2*M_PI;
  }
  else
  {
    if ( pDPhi > 0 )
    {
      fDPhi = pDPhi;
    }
    else
    {
      G4Exception("Error in G4Tubs::G4Tubs - invalid dphi");
    }
  }
	
// Ensure fSphi in 0-2PI or -2PI-0 range if shape crosses 0

  fSPhi = pSPhi;

  if ( fSPhi < 0 )
  {
    fSPhi = 2.0*M_PI - fmod(fabs(fSPhi),2.0*M_PI) ;
  }
  else
  {
    fSPhi = fmod(fSPhi,2.0*M_PI) ;
  }
  if (fSPhi + fDPhi > 2.0*M_PI )
  {
    fSPhi -= 2.0*M_PI ;
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Tubs::~G4Tubs()
{
  ;
}

/////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Tubs::ComputeDimensions(       G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep   )
{
    p->ComputeDimensions(*this,n,pRep) ;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Tubs::CalculateExtent( const EAxis              pAxis,
			        const G4VoxelLimits&     pVoxelLimit,
			        const G4AffineTransform& pTransform,
			              G4double&          pMin, 
                                      G4double&          pMax    ) const
{
  if ( !pTransform.IsRotated() && fDPhi == 2.0*M_PI && fRMin == 0 )
  {
// Special case handling for unrotated solid tubes
// Compute x/y/z mins and maxs fro bounding box respecting limits,
// with early returns if outside limits. Then switch() on pAxis,
// and compute exact x and y limit for x/y case
	    
    G4double xoffset, xMin, xMax ;
    G4double yoffset, yMin, yMax ;
    G4double zoffset, zMin, zMax ;

    G4double diff1, diff2, maxDiff, newMin, newMax ;
    G4double xoff1, xoff2, yoff1, yoff2 ;

    xoffset = pTransform.NetTranslation().x() ;
    xMin = xoffset - fRMax ;
    xMax = xoffset + fRMax ;

    if (pVoxelLimit.IsXLimited())
    {
      if ( xMin > pVoxelLimit.GetMaxXExtent() || 
           xMax < pVoxelLimit.GetMinXExtent()    ) return false;
      else
      {
	if ( xMin < pVoxelLimit.GetMinXExtent() )
	{
	  xMin = pVoxelLimit.GetMinXExtent() ;
	}
	if (xMax > pVoxelLimit.GetMaxXExtent() )
	{
	  xMax = pVoxelLimit.GetMaxXExtent() ;
	}
      }
    }
    yoffset = pTransform.NetTranslation().y() ;
    yMin    = yoffset - fRMax ;
    yMax    = yoffset + fRMax ;

    if ( pVoxelLimit.IsYLimited() )
    {
      if ( yMin > pVoxelLimit.GetMaxYExtent() ||
           yMax < pVoxelLimit.GetMinYExtent()    ) return false ;
      else
      {
	if ( yMin < pVoxelLimit.GetMinYExtent() )
	{
	  yMin = pVoxelLimit.GetMinYExtent() ;
	}
	if ( yMax > pVoxelLimit.GetMaxYExtent() )
	{
				    yMax=pVoxelLimit.GetMaxYExtent();
	}
      }
    }
    zoffset = pTransform.NetTranslation().z() ;
    zMin    = zoffset - fDz ;
    zMax    = zoffset + fDz ;

    if ( pVoxelLimit.IsZLimited() )
    {
      if ( zMin > pVoxelLimit.GetMaxZExtent() ||
           zMax < pVoxelLimit.GetMinZExtent()    ) return false ;
      else
      {
	if ( zMin < pVoxelLimit.GetMinZExtent() )
	{
	  zMin = pVoxelLimit.GetMinZExtent() ;
	}
	if ( zMax > pVoxelLimit.GetMaxZExtent() )
	{
	  zMax=pVoxelLimit.GetMaxZExtent();
	}
      }
    }
    switch ( pAxis )  // Known to cut cylinder
    {
      case kXAxis :
      {
	yoff1 = yoffset - yMin ;
	yoff2 = yMax    - yoffset ;

	if ( yoff1 >= 0 && yoff2 >= 0 ) // Y limits cross max/min x => no change
	{
	  pMin = xMin ;
	  pMax = xMax ;
	}
	else
        {
// Y limits don't cross max/min x => compute max delta x, hence new mins/maxs

	  diff1=sqrt(fRMax*fRMax-yoff1*yoff1);
	  diff2=sqrt(fRMax*fRMax-yoff2*yoff2);
	  maxDiff=(diff1>diff2) ? diff1:diff2;
	  newMin=xoffset-maxDiff;
	  newMax=xoffset+maxDiff;
	  pMin=(newMin<xMin) ? xMin : newMin;
	  pMax=(newMax>xMax) ? xMax : newMax;
	}    
	break;
      }
      case kYAxis :
      {
	xoff1 = xoffset - xMin ;
	xoff2 = xMax - xoffset ;

	if ( xoff1 >= 0 && xoff2 >= 0 ) // X limits cross max/min y => no change
	{
	  pMin = yMin ;
	  pMax = yMax ;
	}
        else
	{
// X limits don't cross max/min y => compute max delta y, hence new mins/maxs

	  diff1   = sqrt(fRMax*fRMax - xoff1*xoff1) ;
	  diff2   = sqrt(fRMax*fRMax - xoff2*xoff2) ;
	  maxDiff = (diff1 > diff2) ? diff1 : diff2 ;
	  newMin  = yoffset - maxDiff ;
	  newMax  = yoffset + maxDiff ;
	  pMin    = (newMin < yMin) ? yMin : newMin ;
	  pMax     =(newMax > yMax) ? yMax : newMax ;
        }
	break ;
      }
      case kZAxis:
      {
	pMin = zMin ;
	pMax = zMax ;
	break ;
      }
    }
    pMin -= kCarTolerance ;
    pMax += kCarTolerance ;
    return true;    
  }
  else // Calculate rotated vertex coordinates
  {
    G4int i, noEntries, noBetweenSections4 ;
    G4bool existsAfterClip = false ;
    G4ThreeVectorList* vertices = CreateRotatedVertices(pTransform) ;

    pMin = +kInfinity ;
    pMax = -kInfinity ;

    noEntries = vertices->entries() ;
    noBetweenSections4 = noEntries - 4 ;
	    
    for (i = 0 ; i < noEntries ; i += 4 )
    {
       ClipCrossSection(vertices,i,pVoxelLimit,pAxis,pMin,pMax) ;
    }
    for (i = 0 ; i < noBetweenSections4 ; i += 4 )
    {
       ClipBetweenSections(vertices,i,pVoxelLimit,pAxis,pMin,pMax) ;
    }
    if (pMin != kInfinity || pMax != -kInfinity )
    {
      existsAfterClip = true ;
      pMin -= kCarTolerance ; // Add 2*tolerance to avoid precision troubles
      pMax += kCarTolerance ;
    }
    else
    {
// Check for case where completely enveloping clipping volume
// If point inside then we are confident that the solid completely
// envelopes the clipping volume. Hence set min/max extents according
// to clipping volume extents along the specified axis.

      G4ThreeVector clipCentre(
      (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
      (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
      (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5 ) ;
		    
      if ( Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside )
      {
	existsAfterClip = true ;
        pMin            = pVoxelLimit.GetMinExtent(pAxis) ;
	pMax            = pVoxelLimit.GetMaxExtent(pAxis) ;
      }
    }
    delete vertices;
    return existsAfterClip;
  }
}

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4Tubs::Inside(const G4ThreeVector& p) const
{
  G4double r2,pPhi,tolRMin,tolRMax;
  EInside in = kOutside ;

  if (fabs(p.z()) <= fDz - kCarTolerance*0.5)
  {
    r2 = p.x()*p.x() + p.y()*p.y() ;

    if (fRMin) tolRMin = fRMin + kRadTolerance*0.5 ;
    else       tolRMin = 0 ;

	       tolRMax = fRMax - kRadTolerance*0.5 ;
	    
    if (r2 >= tolRMin*tolRMin && r2 <= tolRMax*tolRMax)
    {
      if ( fDPhi == 2*M_PI || r2 == 0 )  in = kInside ;
      else
      {
// Try inner tolerant phi boundaries (=>inside)
// if not inside, try outer tolerant phi boundaries

        pPhi = atan2(p.y(),p.x()) ;

	if ( pPhi < 0 ) pPhi += 2*M_PI ; // 0<=pPhi<2pi

	if ( fSPhi >= 0 )
	{
	  if ( pPhi >= fSPhi + kAngTolerance*0.5 &&
               pPhi <= fSPhi + fDPhi-kAngTolerance*0.5 ) in = kInside ;

	  else if ( pPhi >= fSPhi - kAngTolerance*0.5 &&
		    pPhi <= fSPhi + fDPhi + kAngTolerance*0.5 ) in = kSurface ;
	}
	else
	{
	  // if (pPhi < fSPhi + 2*M_PI) pPhi += 2*M_PI ;

	  if ( pPhi >= fSPhi + 2*M_PI + kAngTolerance*0.5 &&
	       pPhi <= fSPhi + fDPhi + 2*M_PI - kAngTolerance*0.5) in = kInside ;

	  else if (pPhi >= fSPhi+2*M_PI-kAngTolerance*0.5 &&
		   pPhi <= fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5) in = kSurface ;
	}			    			    
      }
    }
    else  // Try generous boundaries
    {
      tolRMin = fRMin - kRadTolerance*0.5 ;
      tolRMax = fRMax + kRadTolerance*0.5 ;

      if ( tolRMin < 0 ) tolRMin = 0 ;

      if (r2 >= tolRMin*tolRMin && r2 <= tolRMax*tolRMax)
      {
	if ( fDPhi == 2*M_PI || r2 == 0 ) // Continuous in phi or on z-axis
	{
	    in = kSurface ;
	}
	else // Try outer tolerant phi boundaries only
	{
	  pPhi = atan2(p.y(),p.x()) ;

	  if ( pPhi < 0 ) pPhi += 2*M_PI ; // 0<=pPhi<2pi

	  if ( fSPhi >= 0 )
	  {
	    if (pPhi >= fSPhi-kAngTolerance*0.5 &&
		  pPhi <= fSPhi+fDPhi+kAngTolerance*0.5) in = kSurface ;	
	  }
	  else
	  {
	    //  if (pPhi < fSPhi + 2*M_PI ) pPhi += 2*M_PI ;

	    if (pPhi >= fSPhi+2*M_PI-kAngTolerance*0.5 &&
		pPhi <= fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5) in=kSurface;
	  }
	}
      }
    }
  }
  else if (fabs(p.z()) <= fDz+kCarTolerance*0.5) // Check within tolerant r limits
  {
    r2      = p.x()*p.x() + p.y()*p.y() ;
    tolRMin = fRMin - kRadTolerance*0.5 ;
    tolRMax = fRMax + kRadTolerance*0.5 ;

    if ( tolRMin < 0 ) tolRMin = 0 ;

    if ( r2 >= tolRMin*tolRMin && r2 <= tolRMax*tolRMax )
    {
      if (fDPhi == 2*M_PI || r2 == 0 ) // Continuous in phi or on z-axis
      {
	in = kSurface ;
      }
      else // Try outer tolerant phi boundaries
      {
	pPhi = atan2(p.y(),p.x()) ;

	if ( pPhi < 0 ) pPhi += 2*M_PI ;		// 0<=pPhi<2pi

	if ( fSPhi >= 0 )
	{
	  if (pPhi >= fSPhi-kAngTolerance*0.5 &&
	      pPhi <= fSPhi+fDPhi+kAngTolerance*0.5) in=kSurface; 	
	}
	else
	{
	  //  if ( pPhi < fSPhi + 2*M_PI ) pPhi += 2*M_PI ;

	  if (pPhi>=fSPhi+2*M_PI-kAngTolerance*0.5 &&
		pPhi<=fSPhi+fDPhi+2*M_PI+kAngTolerance*0.5) in = kSurface ;
	}	    
      }
    }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Tubs::SurfaceNormal( const G4ThreeVector& p) const
{
  ENorm side ;
  G4ThreeVector norm ;
  G4double rho, phi ;
  G4double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin ;

  rho = sqrt(p.x()*p.x() + p.y()*p.y()) ;

  distRMin = fabs(rho - fRMin) ;
  distRMax = fabs(rho - fRMax) ;
  distZ    = fabs(fabs(p.z()) - fDz) ;

  if (distRMin < distRMax) // First minimum
  {
    if ( distZ < distRMin )
    {
       distMin = distZ ;
       side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin   ;
    }
  }
  else
  {
    if ( distZ < distRMax )
    {
      distMin = distZ ;
      side    = kNZ   ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax   ;
    }
  }   
  if (fDPhi < 2.0*M_PI  &&  rho ) // Protected against (0,0,z) 
  {
    phi = atan2(p.y(),p.x()) ;

    if ( phi < 0 ) phi += 2*M_PI ;

    if ( fSPhi < 0 )
    {
      distSPhi = fabs(phi - (fSPhi + 2.0*M_PI))*rho ;
    }
    else
    {
	distSPhi=fabs(phi-fSPhi)*rho;
    }
    distEPhi=fabs(phi-fSPhi-fDPhi)*rho;
                                      
    if (distSPhi < distEPhi) // Find new minimum
    {
      if ( distSPhi < distMin )
      {
	   side = kNSPhi ;
      }
    }
    else
    {
      if ( distEPhi < distMin )
      {
	   side=kNEPhi;
      }
    }
  }		
  switch ( side )
  {
    case kNRMin : // Inner radius
    {			                
      norm = G4ThreeVector(-p.x()/rho,-p.y()/rho,0);
      break ;
    }
    case kNRMax : // Outer radius
    {	                
      norm = G4ThreeVector(p.x()/rho,p.y()/rho,0) ;
      break ;
    }
    case kNZ : //    + or - dz
    {			                        
      if ( p.z() > 0 ) norm = G4ThreeVector(0,0,1)  ; 
      else             norm = G4ThreeVector(0,0,-1) ; 
      break ;
    }
    case kNSPhi:
    {
      norm = G4ThreeVector(sin(fSPhi),-cos(fSPhi),0) ;
      break ;
    }
    case kNEPhi:
    {
      norm = G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0) ;
      break;
    }
    default:
    {
      G4Exception("Logic error in G4Tubs::SurfaceNormal") ;
      break ;
    }    
  }                
  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points

G4double G4Tubs::DistanceToIn(const G4ThreeVector& p,
			      const G4ThreeVector& v    ) const
{
  G4double snxt = kInfinity ;  // snxt = default return value

// Precalculated trig for phi intersections - used by r,z intersections to
//                                            check validity

  G4bool seg ;				// true if segmented

  G4double hDPhi, hDPhiOT, hDPhiIT, cosHDPhiOT, cosHDPhiIT ;
					// half dphi + outer tolerance

  G4double cPhi, sinCPhi, cosCPhi ;	// central phi

  G4double tolORMin2, tolIRMax2 ;	// `generous' radii squared

  G4double tolORMax2, tolIRMin2, tolODz, tolIDz ;

// Intersection point variables
  G4double Dist, s, xi, yi, zi, rho2, inum, iden, cosPsi ; 

  G4double t1, t2, t3, b, c, d ;   // Quadratic solver variables 

  G4double Comp ;
  G4double cosSPhi, sinSPhi ;		// Trig for phi start intersect

  G4double ePhi, cosEPhi, sinEPhi ;	// for phi end intersect

// Set phi divided flag and precalcs

  if ( fDPhi < 2.0*M_PI )
  {
     seg=true;
     hDPhi=0.5*fDPhi;		// half delta phi
     cPhi=fSPhi+hDPhi;;
     hDPhiOT=hDPhi+0.5*kAngTolerance;	// outers tol' half delta phi 
     hDPhiIT=hDPhi-0.5*kAngTolerance;
     sinCPhi=sin(cPhi);
     cosCPhi=cos(cPhi);
     cosHDPhiOT=cos(hDPhiOT);
     cosHDPhiIT=cos(hDPhiIT);
  }
  else seg = false ;

// Calculate tolerant rmin and rmax

  if (fRMin > kRadTolerance)
  {
     tolORMin2=(fRMin-0.5*kRadTolerance)*(fRMin-0.5*kRadTolerance);
     tolIRMin2=(fRMin+0.5*kRadTolerance)*(fRMin+0.5*kRadTolerance);
  }
  else
  {
     tolORMin2=0;
     tolIRMin2=0;
  }
  tolORMax2=(fRMax+0.5*kRadTolerance)*(fRMax+0.5*kRadTolerance);
  tolIRMax2=(fRMax-0.5*kRadTolerance)*(fRMax-0.5*kRadTolerance);

// Intersection with Z surfaces

  tolIDz=fDz-kCarTolerance*0.5;
  tolODz=fDz+kCarTolerance*0.5;

  if (fabs(p.z()) >= tolIDz)
  {
     if ( p.z()*v.z() < 0 )	  // at +Z going in -Z or visa versa
     {
	s=(fabs(p.z())-fDz)/fabs(v.z()) ;     // Z intersect distance

        if(s < 0.0) s = 0.0 ;

	xi=p.x()+s*v.x();         	     // Intersection coords
        yi=p.y()+s*v.y();
	rho2=xi*xi+yi*yi;

// Check validity of intersection

	if (tolIRMin2 <= rho2 && rho2 <= tolIRMax2)
	{
	   if (seg && rho2)
	   {
// Psi = angle made with central (average) phi of shape

              inum = xi*cosCPhi+yi*sinCPhi;
              iden = sqrt(rho2);
	      cosPsi = inum/iden;

	      if (cosPsi >= cosHDPhiIT) return s ;
	   }
	   else return s ;
        }
     }
     else
     {
	return snxt ;	// On/outside extent, and heading away
					// -> cannot intersect
     }
  }

// -> Can not intersect z surfaces
//
// Intersection with rmax (possible return) and rmin (must also check phi)
//
// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=R^2
//
// Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
//            t1                t2                t3

  t1=1.0-v.z()*v.z();
  t2=p.x()*v.x()+p.y()*v.y();
  t3=p.x()*p.x()+p.y()*p.y();

  if ( t1 > 0 )				// Check not || to z axis
  {
     b = t2/t1 ;
     c = t3 - fRMax*fRMax ;

     if (t3 >= tolORMax2 && t2<0)   // This also handles the tangent case
     {
// Try outer cylinder intersection
//		      c=(t3-fRMax*fRMax)/t1;

	c /= t1 ;
	d=b*b-c;

	if (d >= 0)	// If real root
	{
	   s = -b-sqrt(d) ;

	   if (s >= 0)	// If 'forwards'
	   {
// Check z intersection
	      zi=p.z()+s*v.z();

	      if (fabs(zi)<=tolODz)
	      {
// Z ok. Check phi intersection if reqd

		  if (!seg) return s ; 
		  else
		  {
		     xi=p.x()+s*v.x();
		     yi=p.y()+s*v.y();
		     cosPsi=(xi*cosCPhi+yi*sinCPhi)/fRMax;

		     if (cosPsi >= cosHDPhiIT) return s ;
		  }
	        }	//  end if fabs(zi)
	     }   	//  end if (s>=0)
	  }	        //  end if (d>=0)
       }	        //  end if (r>=fRMax)

/* *************************************************************************
   // This code is shown here just to make explicit the logic contained
   // in the condition if (t3>tolIRMin2 && t2<0 && fabs(p.z())<=tolIDz),
   // which is placed below.

//	    else if (t3>=tolIRMax2 && t2<0 && fabs(p.z())<=tolIDz)
//	    {                                // Point on Rmax surface
//	       if (!seg)
//	       {
//		  return s = 0 ;   // No Phi cut and move inside
//	       }
//	       else
//	       {
//		  cosPsi=(p.x()*cosCPhi+
//	          p.y()*sinCPhi)/fRMax;
//	          if (cosPsi>=cosHDPhiOT)
//		  { 
//		     return s = 0 ;    // On real Rmax surface and move inside
//		  }
//	       }
//	    }

************************************************************************* */ 

       else 
       {
          // Inside outer radius :
          // check not inside, and heading through tubs (-> 0 to in)

	  if (t3>tolIRMin2 && t2<0 && fabs(p.z())<=tolIDz)
	  {
// Inside both radii, delta r -ve, inside z extent

	     if (seg)
	     {
                inum = p.x()*cosCPhi+p.y()*sinCPhi;
                iden = sqrt(t3);
		cosPsi = inum/iden;

		if (cosPsi >= cosHDPhiIT) return 0 ;
	     }
	     else return 0 ;
	  }
       }	    
       if ( fRMin )    // Try inner cylinder intersection
       {
	  c = (t3 - fRMin*fRMin)/t1 ;
	  d = b*b - c ;

	  if ( d >= 0 )	// If real root
	  {
// Always want 2nd root - we are outside and know rmax Hit was bad
// - If on surface of rmin also need farthest root

	      s = -b + sqrt(d) ;

	      if (s >= -0.5*kCarTolerance)	// check forwards
	      {
// Check z intersection
                 if(s < 0.0) s = 0.0 ;
		 zi=p.z()+s*v.z();

		 if (fabs(zi)<=tolODz)
		 {
// Z ok. Check phi 
		    if (!seg) return s ; 
		    else
		    {
		       xi=p.x()+s*v.x();
		       yi=p.y()+s*v.y();
		       cosPsi=(xi*cosCPhi+yi*sinCPhi)/fRMin;

		       if (cosPsi >= cosHDPhiIT)
		       {
// Good inner radius isect - but earlier phi isect still possible

			  snxt=s;
		       }
		    }
		 }	      //    end if fabs(zi)
	      }	              //    end if (s>=0)
	   }	              //    end if (d>=0)
	}	              //    end if (fRMin)
     }

// Phi segment intersection
//
// o Tolerant of points inside phi planes by up to kCarTolerance*0.5
//
// o NOTE: Large duplication of code between sphi & ephi checks
//         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
//            intersection check <=0 -> >=0
//         -> use some form of loop Construct ?
//
    if ( seg )
    {
// First phi surface (`S'tarting phi)

       sinSPhi=sin(fSPhi);
       cosSPhi=cos(fSPhi);
       Comp=v.x()*sinSPhi-v.y()*cosSPhi;
	                	
       if ( Comp < 0 )  // Component in outwards normal dirn
       {
	  Dist=(p.y()*cosSPhi-p.x()*sinSPhi);

	  if ( Dist < kCarTolerance*0.5 )
	  {
	     s = Dist/Comp ;

	     if (s < snxt)
	     {
		if ( s < 0 ) s = 0 ;
		zi=p.z()+s*v.z();

		if ( fabs(zi) <= tolODz )
		{
		   xi  = p.x() + s*v.x() ;
		   yi  = p.y() + s*v.y() ;
 		   rho2= xi*xi + yi*yi ;

 		   if (   ( rho2 >= tolIRMin2 && rho2 <= tolIRMax2     )

		       || ( rho2 > tolORMin2  && rho2 < tolIRMin2  && 
                            ( v.y()*cosSPhi - v.x()*sinSPhi >  0 ) && 
                            ( v.x()*cosSPhi + v.y()*sinSPhi >= 0 )     )

		       || ( rho2 > tolIRMax2 && rho2 < tolORMax2   && 
                            ( v.y()*cosSPhi - v.x()*sinSPhi > 0 )  && 
                            ( v.x()*cosSPhi + v.y()*sinSPhi < 0 )      )  )
                   {
// z and r intersections good - check intersecting with correct half-plane

		      if ((yi*cosCPhi-xi*sinCPhi) <= 0) snxt = s ;
 		   }    
	        }
	     }
	  }		
       }
	    
// Second phi surface (`E'nding phi)

       ePhi    = fSPhi + fDPhi ;
       sinEPhi = sin(ePhi) ;
       cosEPhi = cos(ePhi) ;
       Comp    = -(v.x()*sinEPhi - v.y()*cosEPhi) ;
				
       if (Comp < 0 )  // Component in outwards normal dirn
       {
	  Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;

	  if ( Dist < kCarTolerance*0.5 )
	  {
	     s = Dist/Comp ;

	     if (s < snxt)
	     {
		if ( s < 0 ) s = 0 ;

		zi=p.z()+s*v.z();

		if ( fabs(zi) <= tolODz )
		{
		   xi   = p.x() + s*v.x() ;
		   yi   = p.y() + s*v.y() ;
 		   rho2 = xi*xi + yi*yi ;

 		   if (   ( rho2 >= tolIRMin2 && rho2 <= tolIRMax2     )

		       || ( rho2 > tolORMin2  && rho2 < tolIRMin2  && 
                            ( v.x()*sinEPhi - v.y()*cosEPhi >  0 ) && 
                            ( v.x()*cosEPhi + v.y()*sinEPhi >= 0 )     )

		       || ( rho2 > tolIRMax2 && rho2 < tolORMax2   && 
                            ( v.x()*sinEPhi - v.y()*cosEPhi > 0 )  && 
                            ( v.x()*cosEPhi + v.y()*sinEPhi < 0 )      )  )
     		   {
// z and r intersections good - check intersecting with correct half-plane

		      if ( (yi*cosCPhi-xi*sinCPhi) >= 0 ) snxt = s ;
 		   }    
		}
	     }
	  }
       }         //  Comp < 0
    }            //  seg != 0
    return snxt;
}
 
//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes 
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4Tubs::DistanceToIn(const G4ThreeVector& p) const
{
  G4double safe, rho, safe1, safe2, safe3 ;
  G4double phiC, cosPhiC, sinPhiC, safePhi, ePhi, cosPsi ;

  rho   = sqrt(p.x()*p.x() + p.y()*p.y()) ;
  safe1 = fRMin - rho ;
  safe2 = rho - fRMax ;
  safe3 = fabs(p.z()) - fDz ;

  if ( safe1 > safe2 ) safe = safe1 ;
  else                 safe = safe2 ;
  if ( safe3 > safe )  safe = safe3 ;

  if (fDPhi < 2.0*M_PI && rho)
  {
    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = cos(phiC) ;
    sinPhiC = sin(phiC) ; 
// Psi=angle from central phi to point
    cosPsi = (p.x()*cosPhiC + p.y()*sinPhiC)/rho ;

    if ( cosPsi < cos(fDPhi*0.5) )
    {
// Point lies outside phi range
      if ( (p.y()*cosPhiC - p.x()*sinPhiC) <= 0 )
      {
	safePhi = fabs(p.x()*sin(fSPhi) - p.y()*cos(fSPhi)) ;
      }
      else
      {
	ePhi    = fSPhi + fDPhi ;
	safePhi = fabs(p.x()*sin(ePhi) - p.y()*cos(ePhi)) ;
      }
      if ( safePhi > safe ) safe = safePhi ;
    }
  }
  if ( safe < 0 ) safe = 0 ;
  return safe ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4Tubs::DistanceToOut( const G4ThreeVector& p,
                                const G4ThreeVector& v,
			        const G4bool calcNorm,
			        G4bool *validNorm,
                                G4ThreeVector *n    ) const
{
  ESide side = kNull , sider, sidephi ;

  G4double snxt, sr, sphi, pdist ;

  G4double deltaR, t1, t2, t3, b, c, d2 ;

// Vars for phi intersection:

  G4double sinSPhi, cosSPhi, ePhi, sinEPhi, cosEPhi ;
  G4double cPhi, sinCPhi, cosCPhi ;
  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi ;

// Z plane intersection

  if (v.z() > 0 )
  {
    pdist = fDz - p.z() ;

    if ( pdist > kCarTolerance*0.5 )
    {
      snxt = pdist/v.z() ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
	*n         = G4ThreeVector(0,0,1) ;
	*validNorm = true ;
      }
      return snxt = 0 ;
    }
  }
  else if ( v.z() < 0 )
  {
    pdist = fDz + p.z() ;

    if ( pdist > kCarTolerance*0.5 )
    {
	 snxt=-pdist/v.z();
	 side=kMZ;
    }
    else
    {
	 if (calcNorm)
	 {
	    *n=G4ThreeVector(0,0,-1);
	    *validNorm=true;
	 }
	 return snxt=0;
    }
  }
  else
  {
      snxt=kInfinity;		// Travel perpendicular to z axis
      side=kNull;
  }

// Radial Intersections
//
// Find intersction with cylinders at rmax/rmin
// Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
//
// Intersects with x^2+y^2=R^2
//
// Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
//
//            t1                t2                    t3

  t1 = 1.0 - v.z()*v.z() ;			// since v normalised
  t2 = p.x()*v.x() + p.y()*v.y() ;
  t3 = p.x()*p.x() + p.y()*p.y() ;

  if ( t1 > 0 )				// Check not parallel
  {
// Calculate sr, r exit distance
	   
    if (t2>=0)
    {
// Delta r not negative => leaving via rmax

      deltaR = t3 - fRMax*fRMax ;

// NOTE: Should use rho-fRMax<-kRadTolerance*0.5 - avoid sqrt for efficiency

      if ( deltaR < -kRadTolerance*fRMax )
      {
	b     = t2/t1 ;
	c     = deltaR/t1 ;
	sr    = -b + sqrt(b*b - c);
        sider = kRMax ;
      }
      else
      {
// On tolerant boundary & heading outwards (or perpendicular to) outer
// radial surface -> leaving immediately

	if (calcNorm) 
	{
	      // if ( p.x() || p.y() )
	      // {
	      //  *n=G4ThreeVector(p.x(),p.y(),0);
	      // }
	      // else
	      // {
	      //  *n=v;
	      // }
	       *n=G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0);
	       *validNorm=true;
	}
	return snxt=0; // Leaving by rmax immediately
      }
    }             
    else  // i.e.  t2 < 0
    {
// Possible rmin intersection

      if ( fRMin )
      {
	deltaR=t3-fRMin*fRMin;
	b=t2/t1;
	c=deltaR/t1;
	d2=b*b-c;

	if ( d2 >= 0 )   // Leaving via rmin
	{
// NOTE: SHould use rho-rmin>kRadTolerance*0.5 - avoid sqrt for efficiency

	  if (deltaR > kRadTolerance*fRMin)
	  {
		  sr = -b-sqrt(d2) ;
		  sider = kRMin ;
	  }
	  else
	  {
	    if ( calcNorm ) *validNorm = false ; // Concave side
	    return snxt=0;
	  }
	}
	else    // No rmin intersect -> must be rmax intersect
	{
	       deltaR=t3-fRMax*fRMax;
	       c=deltaR/t1;
	       sr=-b+sqrt(b*b-c);
	       sider=kRMax;
	}
      }
      else   // No rmin intersect -> must be rmax intersect
      {
	    deltaR=t3-fRMax*fRMax;
	    b=t2/t1;
	    c=deltaR/t1;
	    sr=-b+sqrt(b*b-c);
	    sider=kRMax;
      }
    }
		
// Phi Intersection

	
    if ( fDPhi < 2.0*M_PI )
    {
      sinSPhi=sin(fSPhi);
      cosSPhi=cos(fSPhi);
      ePhi=fSPhi+fDPhi;
      sinEPhi=sin(ePhi);
      cosEPhi=cos(ePhi);
      cPhi=fSPhi+fDPhi*0.5;
      sinCPhi=sin(cPhi);
      cosCPhi=cos(cPhi);


      if ( p.x() || p.y() )  // Check if on z axis (rho not needed later)
      {
// pDist -ve when inside

	pDistS = p.x()*sinSPhi - p.y()*cosSPhi ;
	pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

// Comp -ve when in direction of outwards normal

	 compS   = -sinSPhi*v.x()+cosSPhi*v.y();
	 compE   =  sinEPhi*v.x()-cosEPhi*v.y();
	 sidephi = kNull;

	     //	     if ( pDistS <= 0 && pDistE <= 0 )

         if( ( fDPhi <= pi && ( pDistS <= 0.5*kCarTolerance && 
                                    pDistE <= 0.5*kCarTolerance    ) )   ||
                 ( fDPhi > pi && !( pDistS > 0.5*kCarTolerance && 
                                    pDistE > 0.5*kCarTolerance     ) )  )
	 {
// Inside both phi *full* planes

	   if ( compS < 0 )
	   {
	     sphi=pDistS/compS;
	     xi=p.x()+sphi*v.x();
	     yi=p.y()+sphi*v.y();

// Check intersecting with correct half-plane (if not -> no intersect)

	     if ((yi*cosCPhi-xi*sinCPhi)>=0)
	     {
	       sphi = kInfinity ;
	     }
	     else
	     {
	       sidephi = kSPhi ;

	       if ( pDistS > -kCarTolerance*0.5 )
	       {
			 sphi = 0 ; // Leave by sphi immediately
	       }		
	     }
	   }
	   else sphi = kInfinity ;

	   if ( compE < 0 )
	   {
		   sphi2=pDistE/compE;

// Only check further if < starting phi intersection

		   if (sphi2 > -0.5*kCarTolerance && sphi2 < sphi )
		   {
		      xi=p.x()+sphi2*v.x();
		      yi=p.y()+sphi2*v.y();

// Check intersecting with correct half-plane 

		      if ((yi*cosCPhi-xi*sinCPhi) >= 0)
		      {
// Leaving via ending phi
			 sidephi=kEPhi;

			 if (pDistE<=-kCarTolerance*0.5)
			 {
			    sphi=sphi2;
			 }
			 else 
		         {
			    sphi=0;
			 }
		      }
		   }
		}
	    }
            else sphi = kInfinity ;

	     /* *******************************************

	     else if ( pDistS >= 0 && pDistE >= 0 )
	     {
// Outside both *full* phi planes
                            if (pDistS <= pDistE)
			    {
                              sidephi = kSPhi ;
			    }
                            else
			    {
                              sidephi = kEPhi ;
			    }

		if (fDPhi>M_PI)
		{
		   if (compS<0&&compE<0)
		   {
                      sphi=0;
		   }
		   else 
		   {
                      sphi=kInfinity;
		   }
		}
		else
		{
// if towards both >=0 then once inside (after error) will remain inside

		   if (compS>=0&&compE>=0)
		   {
		      sphi=kInfinity;
		   }
		   else
		   {
		       sphi=0;
		   }
		}
	     }
	     else if ( pDistS > 0 && pDistE < 0 )
	     {
// Outside full starting plane, inside full ending plane

		if (fDPhi>M_PI)
		{
		   if (compE<0)
		   {
		      sphi=pDistE/compE;
		      xi=p.x()+sphi*v.x();
		      yi=p.y()+sphi*v.y();

// Check intersection in correct half-plane (if not -> not leaving phi extent)

		      if ( (yi*cosCPhi-xi*sinCPhi) <= 0 )
		      {
			 sphi=kInfinity;
		      }
		      else   // Leaving via Ending phi
		      {
                         sidephi = kEPhi ;
			 if (pDistE>-kCarTolerance*0.5)
			 {
			    sphi=0;
			 }
		      }
		   }
		   else
		   {
		      sphi=kInfinity;
		   }
		}
		else
	        {
		   if (compS>=0)
		   {
		      if (compE<0)
		      {
			 sphi=pDistE/compE;
			 xi=p.x()+sphi*v.x();
			 yi=p.y()+sphi*v.y();

// Check intersection in correct half-plane (if not -> remain in extent)

			 if ( (yi*cosCPhi-xi*sinCPhi) <= 0 )
			 {
			    sphi=kInfinity;
			 }
			 else  // otherwise leaving via Ending phi
			 {
			    sidephi=kEPhi;
			 }
		      }
		      else 
		      {
                         sphi=kInfinity;
		      }
		   }
		   else  // leaving immediately by starting phi
		   {
		      sidephi=kSPhi;
		      sphi=0;
		   }
		}
	     }
	     else
	     {
// Must be pDistS < 0 && pDistE > 0
// Inside full starting plane, outside full ending plane

		if (fDPhi>M_PI)
		{
		   if (compS<0)
		   {
		      sphi=pDistS/compS;
		      xi=p.x()+sphi*v.x();
		      yi=p.y()+sphi*v.y();

// Check intersection in correct half-plane (if not -> not leaving phi extent)

		      if ((yi*cosCPhi-xi*sinCPhi)>=0)
		      {
			 sphi=kInfinity;
		      }
		      else  // Leaving via Starting phi
		      {
                         sidephi = kSPhi ;   
			 if ( pDistS >- kCarTolerance*0.5 )
			 {
			    sphi=0;
			 }
		      }
		   }
		   else
		   {
		      sphi=kInfinity;
		   }
		}
		else
		{
		   if (compE>=0)
		   {
		      if (compS<0)
		      {
			 sphi=pDistS/compS;
			 xi=p.x()+sphi*v.x();
			 yi=p.y()+sphi*v.y();

// Check intersection in correct half-plane (if not -> remain in extent)

			 if ((yi*cosCPhi-xi*sinCPhi)>=0)
			 {
			    sphi=kInfinity;
			 }
			 else  // otherwise leaving via Starting phi
			 {
			    sidephi=kSPhi;
			 }
		      }
		      else
		      {
			 sphi=kInfinity;
		      }
		   }
		   else  // leaving immediately by ending
		   {
		      sidephi=kEPhi;
		      sphi=0;
		   }
		}
	     }

	     ******************************  */

	  }
	  else
	  {
// On z axis + travel not || to z axis -> if phi of vector direction
// within phi of shape, Step limited by rmax, else Step =0

	     vphi=atan2(v.y(),v.x());

	     if (fSPhi<vphi&&vphi<fSPhi+fDPhi)
	     {
		sphi=kInfinity;
	     }
	     else
	     {
                sidephi = kSPhi ; // arbitrary 
		sphi=0;
	     }
	  }
	  if (sphi<snxt)  // Order intersecttions
	  {
	     snxt=sphi;
	     side=sidephi;
	  }
       }
       if (sr<snxt)  // Order intersections
       {
	  snxt=sr;
	  side=sider;
    }
  }
  if (calcNorm)
  {
    switch(side)
    {
		case kRMax:
					// Note: returned vector not normalised
					// (divide by fRMax for unit vector)
		    xi=p.x()+snxt*v.x();
		    yi=p.y()+snxt*v.y();
		    *n=G4ThreeVector(xi/fRMax,yi/fRMax,0);
		    *validNorm=true;
		    break;

		case kRMin:
		    *validNorm=false;	// Rmin is inconvex
		    break;

		case kSPhi:
		    if (fDPhi<=M_PI)
		    {
		       *n=G4ThreeVector(sin(fSPhi),-cos(fSPhi),0);
		       *validNorm=true;
		    }
		    else
		    {
		       *validNorm=false;
		    }
		    break;

		case kEPhi:
		    if (fDPhi<=M_PI)
		    {
		       *n=G4ThreeVector(-sin(fSPhi+fDPhi),cos(fSPhi+fDPhi),0);
		       *validNorm=true;
		    }
		    else
		    {
		       *validNorm=false;
		    }
		    break;

		case kPZ:
		    *n=G4ThreeVector(0,0,1);
		    *validNorm=true;
		    break;

		case kMZ:
		    *n=G4ThreeVector(0,0,-1);
		    *validNorm=true;
		    break;

		default:
		    G4Exception("Invalid enum in G4Tubs::DistanceToOut");
		    break;
    }
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////////
//
// Calcluate distance (<=actual) to closest surface of shape from inside

G4double G4Tubs::DistanceToOut(const G4ThreeVector& p) const
{
  G4double safe, rho, safeR1, safeR2, safeZ ;
  G4double safePhi, phiC, cosPhiC, sinPhiC, ePhi ;
  rho = sqrt(p.x()*p.x() + p.y()*p.y()) ;

  if ( fRMin )
  {
    safeR1 = rho   - fRMin ;
    safeR2 = fRMax - rho ;
 
    if ( safeR1 < safeR2 ) safe = safeR1 ;
    else                   safe = safeR2 ;
  }
  else safe = fRMax - rho ;
    
  safeZ = fDz - fabs(p.z()) ;

  if ( safeZ < safe ) safe = safeZ ;

// Check if phi divided, Calc distances closest phi plane
  if ( fDPhi < 2.0*M_PI )
  {
// Above/below central phi of Tubs?

    phiC    = fSPhi + fDPhi*0.5 ;
    cosPhiC = cos(phiC) ;
    sinPhiC = sin(phiC) ;

    if ( (p.y()*cosPhiC - p.x()*sinPhiC) <= 0 )
    {
      safePhi = -(p.x()*sin(fSPhi) - p.y()*cos(fSPhi)) ;
    }
    else
    {
      ePhi    = fSPhi + fDPhi ;
      safePhi = (p.x()*sin(ePhi) - p.y()*cos(ePhi)) ;
    }
    if (safePhi < safe) safe = safePhi ;
  }
  if ( safe < 0 ) safe = 0 ;

  return safe ;	
}

/////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

G4ThreeVectorList*
G4Tubs::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
  G4ThreeVectorList* vertices ;
  G4ThreeVector vertex0, vertex1, vertex2, vertex3 ;
  G4double meshAngle, meshRMax, crossAngle, cosCrossAngle, sinCrossAngle, sAngle;
  G4double rMaxX, rMaxY, rMinX, rMinY ;
  G4int crossSection, noCrossSections;

// Compute no of cross-sections necessary to mesh tube
  noCrossSections =G4int(fDPhi/kMeshAngleDefault) + 1 ;

  if ( noCrossSections < kMinMeshSections )
  {
	    noCrossSections = kMinMeshSections;
  }
  else if (noCrossSections>kMaxMeshSections)
  {
	    noCrossSections=kMaxMeshSections;
  }
  meshAngle = fDPhi/(noCrossSections - 1) ;
  meshRMax  = fRMax/cos(meshAngle*0.5) ;

// If complete in phi, set start angle such that mesh will be at fRMax
// on the x axis. Will give better extent calculations when not rotated.

  if (fDPhi == M_PI*2.0 && fSPhi == 0 ) sAngle = -meshAngle*0.5 ;
  else                                  sAngle =  fSPhi ;
    
  vertices=new G4ThreeVectorList(noCrossSections*4);
    
  if ( vertices )
  {
    for (crossSection=0;crossSection<noCrossSections;crossSection++)
    {
// Compute coordinates of cross section at section crossSection
		    crossAngle=sAngle+crossSection*meshAngle;
		    cosCrossAngle=cos(crossAngle);
		    sinCrossAngle=sin(crossAngle);

		    rMaxX=meshRMax*cosCrossAngle;
		    rMaxY=meshRMax*sinCrossAngle;
		    rMinX=fRMin*cosCrossAngle;
		    rMinY=fRMin*sinCrossAngle;
		    vertex0=G4ThreeVector(rMinX,rMinY,-fDz);
		    vertex1=G4ThreeVector(rMaxX,rMaxY,-fDz);
		    vertex2=G4ThreeVector(rMaxX,rMaxY,+fDz);
		    vertex3=G4ThreeVector(rMinX,rMinY,+fDz);

		    vertices->insert(pTransform.TransformPoint(vertex0));
		    vertices->insert(pTransform.TransformPoint(vertex1));
		    vertices->insert(pTransform.TransformPoint(vertex2));
		    vertices->insert(pTransform.TransformPoint(vertex3));
    }
  }
  else
  {
G4Exception("G4Tubs::CreateRotatedVertices Out of memory - Cannot alloc vertices");
  }
  return vertices;
}

///////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation


void G4Tubs::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}

G4Polyhedron* G4Tubs::CreatePolyhedron () const 
{
  return new G4PolyhedronTubs (fRMin, fRMax, fDz, fSPhi, fDPhi);
}

G4NURBS* G4Tubs::CreateNURBS () const 
{
  G4NURBS* pNURBS;
  if (fRMin != 0) 
  {
    if (fDPhi >= 2.0 * M_PI) 
    {
      pNURBS = new G4NURBStube (fRMin, fRMax, fDz);
    }
    else 
    {
      pNURBS = new G4NURBStubesector (fRMin, fRMax, fDz, fSPhi, fSPhi + fDPhi);
    }
  }
  else 
  {
    if (fDPhi >= 2.0 * M_PI) 
    {
      pNURBS = new G4NURBScylinder (fRMax, fDz);
    }
    else 
    {
      const G4double epsilon = 1.e-4; // Cylinder sector not yet available!
      pNURBS = new G4NURBStubesector (epsilon, fRMax, fDz,
				      fSPhi, fSPhi + fDPhi);
    }
  }
  return pNURBS;
}

//
//
/////////////////////////////////// End of G4Tubs.cc ////////////////////////

