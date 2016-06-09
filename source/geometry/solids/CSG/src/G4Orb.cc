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
// $Id: G4Orb.cc,v 1.15 2004/12/02 09:31:28 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// class G4Orb
//
// Implementation for G4Orb class
//
// History:
//
// 30.06.04 V.Grichine - bug fixed in DistanceToIn(p,v) on Rmax surface
// 20.08.03 V.Grichine - created
//
//////////////////////////////////////////////////////////////

#include <assert.h>

#include "G4Orb.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"

// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMax};

// used by normal

enum ENorm {kNRMax};


const G4double G4Orb::fEpsilon = 2.e-11;  // relative tolerance of fRmax

////////////////////////////////////////////////////////////////////////
//
// constructor - check positive radius
//             

G4Orb::G4Orb( const G4String& pName,G4double pRmax )
: G4CSGSolid(pName)
{

  // Check radius

  if (pRmax >= 10*kCarTolerance ) fRmax = pRmax;  
  else
  {
    G4Exception("G4Orb::G4Orb()", "InvalidSetup", FatalException,
                "Invalid radius > 10*kCarTolerance.");
  }
  fRmaxTolerance =  std::max( kRadTolerance, fEpsilon*fRmax);

}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4Orb::~G4Orb()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Orb::ComputeDimensions(       G4VPVParameterisation* p,
                                  const G4int n,
                                  const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*this,n,pRep);
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Orb::CalculateExtent( const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax ) const
{
    // Compute x/y/z mins and maxs for bounding box respecting limits,
    // with early returns if outside limits. Then switch() on pAxis,
    // and compute exact x and y limit for x/y case
      
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;

    G4double diff1,diff2,maxDiff,newMin,newMax;
    G4double xoff1,xoff2,yoff1,yoff2;

    xoffset=pTransform.NetTranslation().x();
    xMin=xoffset-fRmax;
    xMax=xoffset+fRmax;

    if (pVoxelLimit.IsXLimited())
    {
      if ( (xMin>pVoxelLimit.GetMaxXExtent()+kCarTolerance)
        || (xMax<pVoxelLimit.GetMinXExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (xMin<pVoxelLimit.GetMinXExtent())
        {
          xMin=pVoxelLimit.GetMinXExtent();
        }
        if (xMax>pVoxelLimit.GetMaxXExtent())
        {
          xMax=pVoxelLimit.GetMaxXExtent();
        }
      }
    }
    yoffset=pTransform.NetTranslation().y();
    yMin=yoffset-fRmax;
    yMax=yoffset+fRmax;

    if (pVoxelLimit.IsYLimited())
    {
      if ( (yMin>pVoxelLimit.GetMaxYExtent()+kCarTolerance)
        || (yMax<pVoxelLimit.GetMinYExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (yMin<pVoxelLimit.GetMinYExtent())
        {
          yMin=pVoxelLimit.GetMinYExtent();
        }
        if (yMax>pVoxelLimit.GetMaxYExtent())
        {
          yMax=pVoxelLimit.GetMaxYExtent();
        }
      }
    }
    zoffset=pTransform.NetTranslation().z();
    zMin=zoffset-fRmax;
    zMax=zoffset+fRmax;

    if (pVoxelLimit.IsZLimited())
    {
      if ( (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance)
        || (zMax<pVoxelLimit.GetMinZExtent()-kCarTolerance) )
      {
        return false;
      }
      else
      {
        if (zMin<pVoxelLimit.GetMinZExtent())
        {
          zMin=pVoxelLimit.GetMinZExtent();
        }
        if (zMax>pVoxelLimit.GetMaxZExtent())
        {
          zMax=pVoxelLimit.GetMaxZExtent();
        }
      }
    }

    // Known to cut sphere

    switch (pAxis)
    {
      case kXAxis:
        yoff1=yoffset-yMin;
        yoff2=yMax-yoffset;

        if ( yoff1 >= 0 && yoff2 >= 0 )
        {
          // Y limits cross max/min x => no change
          //
          pMin=xMin;
          pMax=xMax;
        }
        else
        {
          // Y limits don't cross max/min x => compute max delta x,
          // hence new mins/maxs
          //
          diff1=std::sqrt(fRmax*fRmax-yoff1*yoff1);
          diff2=std::sqrt(fRmax*fRmax-yoff2*yoff2);
          maxDiff=(diff1>diff2) ? diff1:diff2;
          newMin=xoffset-maxDiff;
          newMax=xoffset+maxDiff;
          pMin=(newMin<xMin) ? xMin : newMin;
          pMax=(newMax>xMax) ? xMax : newMax;
        }
        break;
      case kYAxis:
        xoff1=xoffset-xMin;
        xoff2=xMax-xoffset;
        if (xoff1>=0&&xoff2>=0)
        {
          // X limits cross max/min y => no change
          //
          pMin=yMin;
          pMax=yMax;
        }
        else
        {
          // X limits don't cross max/min y => compute max delta y,
          // hence new mins/maxs
          //
          diff1=std::sqrt(fRmax*fRmax-xoff1*xoff1);
          diff2=std::sqrt(fRmax*fRmax-xoff2*xoff2);
          maxDiff=(diff1>diff2) ? diff1:diff2;
          newMin=yoffset-maxDiff;
          newMax=yoffset+maxDiff;
          pMin=(newMin<yMin) ? yMin : newMin;
          pMax=(newMax>yMax) ? yMax : newMax;
        }
        break;
      case kZAxis:
        pMin=zMin;
        pMax=zMax;
        break;
      default:
        break;
    }
    pMin -= fRmaxTolerance;
    pMax += fRmaxTolerance;

    return true;  
  
}

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
// Split into radius checks
// 

EInside G4Orb::Inside( const G4ThreeVector& p ) const
{
  G4double rad2,tolRMax;
  EInside in;


  rad2 = p.x()*p.x()+p.y()*p.y()+p.z()*p.z() ;

  // G4double rad = std::sqrt(rad2);
  // Check radial surface
  // sets `in'
  
  tolRMax = fRmax - fRmaxTolerance*0.5 ;
    
  if ( rad2 <= tolRMax*tolRMax )  in = kInside ;
  else
  {
    tolRMax = fRmax + fRmaxTolerance*0.5 ;       
    if ( rad2 <= tolRMax*tolRMax ) in = kSurface ;
    else                           in = kOutside ;
  }
  return in;
}

/////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4Orb::SurfaceNormal( const G4ThreeVector& p ) const
{
  ENorm side = kNRMax;
  G4ThreeVector norm;
  G4double rad = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());

  switch (side)
  {
    case kNRMax: 
      norm = G4ThreeVector(p.x()/rad,p.y()/rad,p.z()/rad);
      break;
   default:
      DumpInfo();
      G4Exception("G4Orb::SurfaceNormal()", "Notification", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break;    
  } 

  return norm;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// -> If point is outside outer radius, compute intersection with rmax
//        - if no intersection return
//        - if  valid phi,theta return intersection Dist

G4double G4Orb::DistanceToIn( const G4ThreeVector& p,
                              const G4ThreeVector& v  ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value

  G4double rad2, pDotV3d, tolORMax2, tolIRMax2 ;
  G4double c, d2, s = kInfinity ;

  // General Precalcs

  rad2    = p.x()*p.x() + p.y()*p.y() + p.z()*p.z() ;
  pDotV3d = p.x()*v.x() + p.y()*v.y() + p.z()*v.z() ;

  // Radial Precalcs

  tolORMax2 = (fRmax+fRmaxTolerance*0.5)*(fRmax+fRmaxTolerance*0.5) ;
  tolIRMax2 = (fRmax-fRmaxTolerance*0.5)*(fRmax-fRmaxTolerance*0.5) ;

  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer G4Orb heading through solid (-> 0)
  // - No intersect -> no intersection with G4Orb
  //
  // Shell eqn: x^2+y^2+z^2 = RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2s(pDotV3d)       +s^2                =R^2
  //
  // => s=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  c = rad2 - fRmax*fRmax ;

  if ( c > fRmaxTolerance*fRmax )
  {
    // If outside tolerant boundary of outer G4Orb
    // [ should be std::sqrt(rad2) - fRmax > fRmaxTolerance*0.5 ]

    d2 = pDotV3d*pDotV3d - c ;

    if ( d2 >= 0 )
    {
      s = -pDotV3d - std::sqrt(d2) ;

      if (s >= 0 ) return snxt = s;
           
    }
    else    // No intersection with G4Orb
    {
      return snxt = kInfinity;
    }
  }
  else
  {
    if ( c > -fRmaxTolerance*fRmax )  // on surface  
    {
      d2 = pDotV3d*pDotV3d - c ;             
      //  if ( pDotV3d >= 0 ) return snxt = kInfinity;
      if ( d2 < fRmaxTolerance*fRmax || pDotV3d >= 0 ) return snxt = kInfinity;
      else                return snxt = 0.;
    }
    else // inside ???
    {
      G4Exception("G4Orb::DistanceToIn(p,v)", "Notification",
                  JustWarning, "Point p is inside !?");
    }
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to radial plane
// - Return 0 if point inside

G4double G4Orb::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe=0.0, rad  = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
                 safe = rad - fRmax;
  if( safe < 0 ) safe = 0. ;
  return safe;
}

/////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// 

G4double G4Orb::DistanceToOut( const G4ThreeVector& p,
                               const G4ThreeVector& v,
                               const G4bool calcNorm,
                                     G4bool *validNorm,
                                     G4ThreeVector *n   ) const
{
  G4double snxt = kInfinity;     // ??? snxt is default return value
  ESide    side = kNull;
  
  G4double rad2,pDotV3d; 
  G4double xi,yi,zi;      // Intersection point
  G4double c,d2;
                 
  rad2    = p.x()*p.x() + p.y()*p.y() + p.z()*p.z();
  pDotV3d = p.x()*v.x() + p.y()*v.y() + p.z()*v.z();
    
  // Radial Intersection from G4Orb::DistanceToIn
  //
  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer G4Orb heading through solid (-> 0)
  // - No intersect -> no intersection with G4Orb
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2s(pxvx+pyvy+pzvz)+s^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2s(pDotV3d)       +s^2                =R^2
  //
  // => s=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))
  
  const G4double  Rmax_plus = fRmax + fRmaxTolerance*0.5;

  if( rad2 <= Rmax_plus*Rmax_plus )
  {
    c = rad2-fRmax*fRmax ;

    if ( c < fRmaxTolerance*fRmax) 
    {
      // Within tolerant Outer radius 
      // 
      // The test is
      //     rad  - fRmax < 0.5*fRmaxTolerance
      // =>  rad  < fRmax + 0.5*kRadTol
      // =>  rad2 < (fRmax + 0.5*kRadTol)^2
      // =>  rad2 < fRmax^2 + 2.*0.5*fRmax*kRadTol + 0.25*kRadTol*kRadTol
      // =>  rad2 - fRmax^2    <~    fRmax*kRadTol 

      d2 = pDotV3d*pDotV3d - c;

      if( ( c > -fRmaxTolerance*fRmax) &&         // on tolerant surface
	  ( ( pDotV3d >= 0 )   || ( d2 < 0 )) )   // leaving outside from Rmax 
                                                  // not re-entering
      {
        if(calcNorm)
        {
          *validNorm = true ;
          *n         = G4ThreeVector(p.x()/fRmax,p.y()/fRmax,p.z()/fRmax) ;
        }
        return snxt = 0;
      }
      else 
      {
        snxt = -pDotV3d + std::sqrt(d2);    // second root since inside Rmax
        side = kRMax ; 
      }
    }
  }
  else // p is outside ???
  {
    G4cout.precision(16);
    G4cout << G4endl;
    DumpInfo();
    G4cout << "Position:"  << G4endl << G4endl;
    G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
    G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
    G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
    G4cout << "Rp = "<< std::sqrt( p.x()*p.x()+p.y()*p.y()+p.z()*p.z() )/mm << " mm" 
           << G4endl << G4endl;
    G4cout << "Direction:" << G4endl << G4endl;
    G4cout << "v.x() = "   << v.x() << G4endl;
    G4cout << "v.y() = "   << v.y() << G4endl;
    G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
    G4cout << "Proposed distance :" << G4endl << G4endl;
    G4cout << "snxt = "    << snxt/mm << " mm" << G4endl << G4endl;
    G4Exception("G4Orb::DistanceToOut(p,v,..)", "Notification",
                JustWarning, "Logic error: snxt = kInfinity ???");
  }
  if (calcNorm)    // Output switch operator
  {
    switch( side )
    {
      case kRMax:
        xi=p.x()+snxt*v.x();
        yi=p.y()+snxt*v.y();
        zi=p.z()+snxt*v.z();
        *n=G4ThreeVector(xi/fRmax,yi/fRmax,zi/fRmax);
        *validNorm=true;
        break;
      default:
        G4cout.precision(16);
        G4cout << G4endl;
        DumpInfo();
        G4cout << "Position:"  << G4endl << G4endl;
        G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
        G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
        G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
        G4cout << "Direction:" << G4endl << G4endl;
        G4cout << "v.x() = "   << v.x() << G4endl;
        G4cout << "v.y() = "   << v.y() << G4endl;
        G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
        G4cout << "Proposed distance :" << G4endl << G4endl;
        G4cout << "snxt = "    << snxt/mm << " mm" << G4endl << G4endl;
        G4Exception("G4Orb::DistanceToOut(p,v,..)","Notification",JustWarning,
                    "Undefined side for valid surface normal to solid.");
        break;
    }
  }
  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance (<=actual) to closest surface of shape from inside

G4double G4Orb::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0,rad = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4Exception("G4Orb::DistanceToOut(p)", "Notification", JustWarning, 
                 "Point p is outside !?" );
  }
#endif

  safe = fRmax - rad;
  if ( safe < 0. ) safe = 0.;
  return safe;
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType

G4GeometryType G4Orb::GetEntityType() const
{
  return G4String("G4Orb");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Orb::StreamInfo( std::ostream& os ) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Orb\n"
     << " Parameters: \n"

     << "    outer radius: " << fRmax/mm << " mm \n"
     << "-----------------------------------------------------------\n";

  return os;
}

/////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Orb::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddThis (*this);
}

G4Polyhedron* G4Orb::CreatePolyhedron () const
{
  return new G4PolyhedronSphere (0., fRmax, 0., 2*pi, 0., pi);
}

G4NURBS* G4Orb::CreateNURBS () const
{
  return new G4NURBSbox (fRmax, fRmax, fRmax);       // Box for now!!!
}

// End of G4Orb.cc 
//
/////////////////////////////////////////////////////////////////////
