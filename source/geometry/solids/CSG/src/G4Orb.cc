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
// $Id: G4Orb.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// class G4Orb
//
// Implementation for G4Orb class
//
// History:
//
// 27.10.16 E.Tcherniaev - reimplemented CalculateExtent()
// 05.04.12 M.Kelsey   - GetPointOnSurface() throw flat in cos(theta)
// 30.06.04 V.Grichine - bug fixed in DistanceToIn(p,v) on Rmax surface
// 20.08.03 V.Grichine - created
//
//////////////////////////////////////////////////////////////

#include "G4Orb.hh"

#if !defined(G4GEOM_USE_UORB)

#include "G4TwoVector.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "Randomize.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"

using namespace CLHEP;

// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull,kRMax};

// used by normal

enum ENorm {kNRMax};


////////////////////////////////////////////////////////////////////////
//
// constructor - check positive radius
//             

G4Orb::G4Orb( const G4String& pName, G4double pRmax )
: G4CSGSolid(pName), fRmax(pRmax)
{

  const G4double fEpsilon = 2.e-11;  // relative tolerance of fRmax

  G4double kRadTolerance
    = G4GeometryTolerance::GetInstance()->GetRadialTolerance();

  // Check radius
  //
  if ( pRmax < 10*kCarTolerance )
  {
    G4Exception("G4Orb::G4Orb()", "GeomSolids0002", FatalException,
                "Invalid radius < 10*kCarTolerance.");
  }
  fRmaxTolerance =  std::max( kRadTolerance, fEpsilon*fRmax);

}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Orb::G4Orb( __void__& a )
  : G4CSGSolid(a), fRmax(0.), fRmaxTolerance(0.)
{
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

G4Orb::~G4Orb()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4Orb::G4Orb(const G4Orb& rhs)
  : G4CSGSolid(rhs), fRmax(rhs.fRmax), fRmaxTolerance(rhs.fRmaxTolerance)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4Orb& G4Orb::operator = (const G4Orb& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   fRmax = rhs.fRmax;
   fRmaxTolerance = rhs.fRmaxTolerance;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Orb::ComputeDimensions(       G4VPVParameterisation* p,
                               const G4int n,
                               const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Orb::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  G4double radius = GetRadius();
  pMin.set(-radius,-radius,-radius);
  pMax.set( radius, radius, radius);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4Orb::BoundingLimits()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Orb::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Check bounding box
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Find bounding envelope and calculate extent
  //
  static const G4int NTHETA = 8;  // number of steps along Theta
  static const G4int NPHI   = 16; // number of steps along Phi
  static const G4double sinHalfTheta = std::sin(halfpi/NTHETA);
  static const G4double cosHalfTheta = std::cos(halfpi/NTHETA);
  static const G4double sinHalfPhi   = std::sin(pi/NPHI);
  static const G4double cosHalfPhi   = std::cos(pi/NPHI);
  static const G4double sinStepTheta = 2.*sinHalfTheta*cosHalfTheta;
  static const G4double cosStepTheta = 1. - 2.*sinHalfTheta*sinHalfTheta;
  static const G4double sinStepPhi   = 2.*sinHalfPhi*cosHalfPhi;
  static const G4double cosStepPhi   = 1. - 2.*sinHalfPhi*sinHalfPhi;

  G4double radius = GetRadius();
  G4double rtheta = radius/cosHalfTheta;
  G4double rphi   = rtheta/cosHalfPhi;

  // set reference circle
  G4TwoVector xy[NPHI];
  G4double sinCurPhi = sinHalfPhi;
  G4double cosCurPhi = cosHalfPhi;
  for (G4int k=0; k<NPHI; ++k)
  {
    xy[k].set(cosCurPhi,sinCurPhi);
    G4double sinTmpPhi = sinCurPhi;
    sinCurPhi = sinCurPhi*cosStepPhi + cosCurPhi*sinStepPhi;
    cosCurPhi = cosCurPhi*cosStepPhi - sinTmpPhi*sinStepPhi;
  }
  
  // set bounding circles
  G4ThreeVectorList circles[NTHETA];
  for (G4int i=0; i<NTHETA; ++i) circles[i].resize(NPHI);

  G4double sinCurTheta = sinHalfTheta;
  G4double cosCurTheta = cosHalfTheta;
  for (G4int i=0; i<NTHETA; ++i)
  {
    G4double z = rtheta*cosCurTheta;
    G4double rho = rphi*sinCurTheta;
    for (G4int k=0; k<NPHI; ++k)
    {
      circles[i][k].set(rho*xy[k].x(),rho*xy[k].y(),z);
    }
    G4double sinTmpTheta = sinCurTheta;
    sinCurTheta = sinCurTheta*cosStepTheta + cosCurTheta*sinStepTheta;
    cosCurTheta = cosCurTheta*cosStepTheta - sinTmpTheta*sinStepTheta;
  }

  // set envelope and calculate extent
  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(NTHETA);
  for (G4int i=0; i<NTHETA; ++i) polygons[i] = &circles[i];

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
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


  rad2 = p.x()*p.x()+p.y()*p.y()+p.z()*p.z();

  G4double radius = std::sqrt(rad2);

  // G4double radius = std::sqrt(rad2);
  // Check radial surface
  // sets `in'
  
  tolRMax = fRmax - fRmaxTolerance*0.5;
    
  if ( radius <= tolRMax )  { in = kInside; }
  else
  {
    tolRMax = fRmax + fRmaxTolerance*0.5;       
    if ( radius <= tolRMax )  { in = kSurface; }
    else                   { in = kOutside; }
  }
  return in;
}

/////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p

G4ThreeVector G4Orb::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4double radius = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
  return G4ThreeVector(p.x()/radius,p.y()/radius,p.z()/radius);
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
  G4double snxt = kInfinity;      // snxt = default return value

  G4double radius, pDotV3d; // , tolORMax2, tolIRMax2;
  G4double c, d2, sd = kInfinity;

  const G4double dRmax = 100.*fRmax;

  // General Precalcs

  radius  = std::sqrt(p.x()*p.x() + p.y()*p.y() + p.z()*p.z());
  pDotV3d = p.x()*v.x() + p.y()*v.y() + p.z()*v.z();

  // Radial Precalcs

  // tolORMax2 = (fRmax+fRmaxTolerance*0.5)*(fRmax+fRmaxTolerance*0.5);
  // tolIRMax2 = (fRmax-fRmaxTolerance*0.5)*(fRmax-fRmaxTolerance*0.5);

  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer G4Orb heading through solid (-> 0)
  // - No intersect -> no intersection with G4Orb
  //
  // Shell eqn: x^2+y^2+z^2 = RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2sd(pxvx+pyvy+pzvz)+sd^2(vx^2+vy^2+vz^2)=R^2
  // =>      rad2        +2sd(pDotV3d)      +sd^2                =R^2
  //
  // => sd=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  c = (radius - fRmax)*(radius + fRmax);

  if( radius > fRmax-fRmaxTolerance*0.5 ) // not inside in terms of Inside(p)
  {
    if ( c > fRmaxTolerance*fRmax )
    {
      // If outside tolerant boundary of outer G4Orb in terms of c
      // [ should be std::sqrt(rad2) - fRmax > fRmaxTolerance*0.5 ]

      d2 = pDotV3d*pDotV3d - c;

      if ( d2 >= 0 )
      {
        sd = -pDotV3d - std::sqrt(d2);
        if ( sd >= 0 )
        {
          if ( sd > dRmax ) // Avoid rounding errors due to precision issues seen on
          {                 // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd - std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          } 
          return snxt = sd;
        }
      }
      else    // No intersection with G4Orb
      {
        return snxt = kInfinity;
      }
    }
    else // not outside in terms of c
    {
      if ( c > -fRmaxTolerance*fRmax )  // on surface  
      {
        d2 = pDotV3d*pDotV3d - c;             
        if ( (d2 < fRmaxTolerance*fRmax) || (pDotV3d >= 0) )
        {
          return snxt = kInfinity;
        }
        else
        {
          return snxt = 0.;
        }
      }
    }
  }
#ifdef G4CSGDEBUG
  else // inside ???
  {
      G4Exception("G4Orb::DistanceToIn(p,v)", "GeomSolids1002",
                  JustWarning, "Point p is inside !?");
  }
#endif

  return snxt;
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to radial plane
// - Return 0 if point inside

G4double G4Orb::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safe = 0.0,
           radius  = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());
  safe = radius - fRmax;
  if( safe < 0 ) { safe = 0.; }
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
  G4double radius = std::sqrt(rad2);

  if ( radius <= Rmax_plus )
  {
    c = (radius - fRmax)*(radius + fRmax);

    if ( c < fRmaxTolerance*fRmax ) 
    {
      // Within tolerant Outer radius 
      // 
      // The test is
      //     radius  - fRmax < 0.5*fRmaxTolerance
      // =>  radius  < fRmax + 0.5*kRadTol
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
          *validNorm = true;
          *n         = G4ThreeVector(p.x()/fRmax,p.y()/fRmax,p.z()/fRmax);
        }
        return snxt = 0;
      }
      else 
      {
        snxt = -pDotV3d + std::sqrt(d2);    // second root since inside Rmax
        side = kRMax; 
      }
    }
  }
  else // p is outside ???
  {
    G4cout << G4endl;
    DumpInfo();
    std::ostringstream message;
    G4int oldprc = message.precision(16);
    message << "Logic error: snxt = kInfinity ???" << G4endl
            << "Position:"  << G4endl << G4endl
            << "p.x() = "   << p.x()/mm << " mm" << G4endl
            << "p.y() = "   << p.y()/mm << " mm" << G4endl
            << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
            << "Rp = "<< std::sqrt( p.x()*p.x()+p.y()*p.y()+p.z()*p.z() )/mm
            << " mm" << G4endl << G4endl
            << "Direction:" << G4endl << G4endl
            << "v.x() = "   << v.x() << G4endl
            << "v.y() = "   << v.y() << G4endl
            << "v.z() = "   << v.z() << G4endl << G4endl
            << "Proposed distance :" << G4endl << G4endl
            << "snxt = "    << snxt/mm << " mm" << G4endl;
    message.precision(oldprc);
    G4Exception("G4Orb::DistanceToOut(p,v,..)", "GeomSolids1002",
                JustWarning, message);
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
        G4cout << G4endl;
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl;
        message.precision(oldprc);
        G4Exception("G4Orb::DistanceToOut(p,v,..)","GeomSolids1002",
                    JustWarning, message);
        break;
    }
  }
  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4Orb::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safe=0.0,radius = std::sqrt(p.x()*p.x()+p.y()*p.y()+p.z()*p.z());

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4int oldprc = G4cout.precision(16);
     G4cout << G4endl;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
     G4cout.precision(oldprc);
     G4Exception("G4Orb::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, "Point p is outside !?" );
  }
#endif

  safe = fRmax - radius;
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
// Make a clone of the object
//
G4VSolid* G4Orb::Clone() const
{
  return new G4Orb(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Orb::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Orb\n"
     << " Parameters: \n"

     << "    outer radius: " << fRmax/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4Orb::GetPointOnSurface() const
{
  //  generate a random number from zero to 2pi...
  //
  G4double phi      = G4RandFlat::shoot(0.,2.*pi);
  G4double cosphi   = std::cos(phi);
  G4double sinphi   = std::sin(phi);

  // generate a random point uniform in area
  G4double costheta = G4RandFlat::shoot(-1.,1.);
  G4double sintheta = std::sqrt(1.-sqr(costheta));
  
  return G4ThreeVector (fRmax*sintheta*cosphi,
                        fRmax*sintheta*sinphi, fRmax*costheta); 
}

////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Orb::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this);
}

G4Polyhedron* G4Orb::CreatePolyhedron () const
{
  return new G4PolyhedronSphere (0., fRmax, 0., 2*pi, 0., pi);
}

#endif
