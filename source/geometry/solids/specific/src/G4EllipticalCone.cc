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
// Implementation of G4EllipticalCone class
//
// This code implements an Elliptical Cone given explicitly by the
// equation:
//   x^2/a^2 + y^2/b^2 = (z-h)^2
// and specified by the parameters (a,b,h) and a cut parallel to the
// xy plane above z = 0.
//
// Author: Dionysios Anninos
// Revised: Evgueni Tcherniaev
// --------------------------------------------------------------------

#if !(defined(G4GEOM_USE_UELLIPTICALCONE) && defined(G4GEOM_USE_SYS_USOLIDS))

#include "globals.hh"

#include "G4EllipticalCone.hh"

#include "G4RandomTools.hh"
#include "G4GeomTools.hh"
#include "G4ClippablePolygon.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"

#include "meshdefs.hh"

#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters

G4EllipticalCone::G4EllipticalCone(const G4String& pName,
                                         G4double  pxSemiAxis,
                                         G4double  pySemiAxis,
                                         G4double  pzMax,
                                         G4double  pzTopCut)
  : G4VSolid(pName), zTopCut(0.)
{
  halfCarTol = 0.5*kCarTolerance;

  // Check Semi-Axis & Z-cut
  //
  if ( (pxSemiAxis <= 0.) || (pySemiAxis <= 0.) || (pzMax <= 0.) )
  {
     std::ostringstream message;
     message << "Invalid semi-axis or height for solid: " << GetName()
             << "\n   X semi-axis, Y semi-axis, height = "
             << pxSemiAxis << ", " << pySemiAxis << ", " << pzMax;
     G4Exception("G4EllipticalCone::G4EllipticalCone()", "GeomSolids0002",
                 FatalErrorInArgument, message);
   }

  if ( pzTopCut <= 0 )
  {
     std::ostringstream message;
     message << "Invalid z-coordinate for cutting plane for solid: " << GetName()
             << "\n   Z top cut = " << pzTopCut;
     G4Exception("G4EllipticalCone::G4EllipticalCone()", "GeomSolids0002",
                 FatalErrorInArgument, message);
  }

  SetSemiAxis( pxSemiAxis, pySemiAxis, pzMax );
  SetZCut(pzTopCut);
}

/////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4EllipticalCone::G4EllipticalCone( __void__& a )
  : G4VSolid(a), halfCarTol(0.),
    xSemiAxis(0.), ySemiAxis(0.), zheight(0.), zTopCut(0.),
    cosAxisMin(0.), invXX(0.), invYY(0.)
{
}

/////////////////////////////////////////////////////////////////////////
//
// Destructor

G4EllipticalCone::~G4EllipticalCone()
{
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

/////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4EllipticalCone::G4EllipticalCone(const G4EllipticalCone& rhs)
  : G4VSolid(rhs), halfCarTol(rhs.halfCarTol),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    xSemiAxis(rhs.xSemiAxis), ySemiAxis(rhs.ySemiAxis),
    zheight(rhs.zheight), zTopCut(rhs.zTopCut),
    cosAxisMin(rhs.cosAxisMin), invXX(rhs.invXX), invYY(rhs.invYY)
{
}

/////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4EllipticalCone& G4EllipticalCone::operator = (const G4EllipticalCone& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   halfCarTol = rhs.halfCarTol;
   fCubicVolume = rhs.fCubicVolume; fSurfaceArea = rhs.fSurfaceArea;
   xSemiAxis = rhs.xSemiAxis; ySemiAxis = rhs.ySemiAxis;
   zheight = rhs.zheight; zTopCut = rhs.zTopCut;
   cosAxisMin = rhs.cosAxisMin; invXX = rhs.invXX; invYY = rhs.invYY;

   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

/////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4EllipticalCone::BoundingLimits(G4ThreeVector& pMin,
                                      G4ThreeVector& pMax) const
{
  G4double zcut   = GetZTopCut();
  G4double height = GetZMax(); 
  G4double xmax   = GetSemiAxisX()*(height+zcut);
  G4double ymax   = GetSemiAxisY()*(height+zcut);
  pMin.set(-xmax,-ymax,-zcut);
  pMax.set( xmax, ymax, zcut);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4EllipticalCone::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4EllipticalCone::CalculateExtent(const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin,bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  static const G4int NSTEPS = 48; // number of steps for whole circle
  static const G4double ang = twopi/NSTEPS;
  static const G4double sinHalf = std::sin(0.5*ang);
  static const G4double cosHalf = std::cos(0.5*ang);
  static const G4double sinStep = 2.*sinHalf*cosHalf;
  static const G4double cosStep = 1. - 2.*sinHalf*sinHalf;
  G4double zcut   = bmax.z();
  G4double height = GetZMax(); 
  G4double sxmin  = GetSemiAxisX()*(height-zcut)/cosHalf;
  G4double symin  = GetSemiAxisY()*(height-zcut)/cosHalf;
  G4double sxmax  = bmax.x()/cosHalf;
  G4double symax  = bmax.y()/cosHalf;

  G4double sinCur = sinHalf;
  G4double cosCur = cosHalf;
  G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
  for (G4int k=0; k<NSTEPS; ++k)
  {
    baseA[k].set(sxmax*cosCur,symax*sinCur,-zcut);
    baseB[k].set(sxmin*cosCur,symin*sinCur, zcut);
    
    G4double sinTmp = sinCur;
    sinCur = sinCur*cosStep + cosCur*sinStep;
    cosCur = cosCur*cosStep - sinTmp*sinStep;
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;
  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

/////////////////////////////////////////////////////////////////////////
//
// Determine where is point: inside, outside or on surface

EInside G4EllipticalCone::Inside(const G4ThreeVector& p) const
{
  G4double hp = std::sqrt(p.x()*p.x()*invXX + p.y()*p.y()*invYY) + p.z();
  G4double ds = (hp - zheight)*cosAxisMin;
  G4double dz = std::abs(p.z()) - zTopCut;
  G4double dist = std::max(ds,dz);

  if (dist > halfCarTol) return kOutside;
  return (dist > -halfCarTol) ? kSurface : kInside;
}

/////////////////////////////////////////////////////////////////////////
//
// Return unit normal at surface closest to p

G4ThreeVector G4EllipticalCone::SurfaceNormal( const G4ThreeVector& p) const
{
  G4ThreeVector norm(0,0,0);
  G4int nsurf = 0;  // number of surfaces where p is placed

  G4double hp = std::sqrt(p.x()*p.x()*invXX + p.y()*p.y()*invYY) + p.z();
  G4double ds = (hp - zheight)*cosAxisMin;
  if (std::abs(ds) <= halfCarTol)
  {
    norm = G4ThreeVector(p.x()*invXX, p.y()*invYY, hp - p.z());
    G4double mag = norm.mag();
    if (mag == 0) return G4ThreeVector(0,0,1); // apex
    norm *= (1/mag);
    ++nsurf;
  }
  G4double dz = std::abs(p.z()) - zTopCut;
  if (std::abs(dz) <= halfCarTol)
  {
    norm += G4ThreeVector(0., 0.,(p.z() < 0) ? -1. : 1.);
    ++nsurf;
  }

  if      (nsurf == 1) return norm;
  else if (nsurf >  1) return norm.unit(); // elliptic edge
  else
  {
    // Point is not on the surface
    //
#ifdef G4CSGDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << G4endl;
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4EllipticalCone::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

/////////////////////////////////////////////////////////////////////////
//
// Find surface nearest to point and return corresponding normal.
// The algorithm is similar to the algorithm used in Inside().
// This method normally should not be called.

G4ThreeVector
G4EllipticalCone::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  G4double hp = std::sqrt(p.x()*p.x()*invXX + p.y()*p.y()*invYY) + p.z();
  G4double ds = (hp - zheight)*cosAxisMin;
  G4double dz = std::abs(p.z()) - zTopCut;
  if (ds > dz && std::abs(hp - p.z()) > halfCarTol)
    return G4ThreeVector(p.x()*invXX, p.y()*invYY, hp - p.z()).unit();
  else
    return G4ThreeVector(0., 0.,(p.z() < 0) ? -1. : 1.);
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// return kInfinity if no intersection, or intersection distance <= tolerance

G4double G4EllipticalCone::DistanceToIn( const G4ThreeVector& p,
                                         const G4ThreeVector& v  ) const
{
  G4double distMin = kInfinity;

  // code from EllipticalTube

  G4double sigz = p.z()+zTopCut;

  //
  // Check z = -dz planer surface
  //

  if (sigz < halfCarTol)
  {
    //
    // We are "behind" the shape in z, and so can
    // potentially hit the rear face. Correct direction?
    //
    if (v.z() <= 0)
    {
      //
      // As long as we are far enough away, we know we
      // can't intersect
      //
      if (sigz < 0) return kInfinity;
      
      //
      // Otherwise, we don't intersect unless we are
      // on the surface of the ellipse
      //

      if ( sqr(p.x()/( xSemiAxis - halfCarTol ))
         + sqr(p.y()/( ySemiAxis - halfCarTol )) <= sqr( zheight + zTopCut ) )
        return kInfinity;

    }
    else
    {
      //
      // How far?
      //
      G4double q = -sigz/v.z();
      
      //
      // Where does that place us?
      //
      G4double xi = p.x() + q*v.x(),
               yi = p.y() + q*v.y();
      
      //
      // Is this on the surface (within ellipse)?
      //
      if ( sqr(xi/xSemiAxis) + sqr(yi/ySemiAxis) <= sqr( zheight + zTopCut ) )
      {
        //
        // Yup. Return q, unless we are on the surface
        //
        return (sigz < -halfCarTol) ? q : 0;
      }
      else if (xi/(xSemiAxis*xSemiAxis)*v.x()
             + yi/(ySemiAxis*ySemiAxis)*v.y() >= 0)
      {
        //
        // Else, if we are traveling outwards, we know
        // we must miss
        //
        //        return kInfinity;
      }
    }
  }

  //
  // Check z = +dz planer surface
  //
  sigz = p.z() - zTopCut;
  
  if (sigz > -halfCarTol)
  {
    if (v.z() >= 0)
    {

      if (sigz > 0) return kInfinity;

      if ( sqr(p.x()/( xSemiAxis - halfCarTol ))
         + sqr(p.y()/( ySemiAxis - halfCarTol )) <= sqr( zheight-zTopCut ) )
        return kInfinity;

    }
    else {
      G4double q = -sigz/v.z();

      G4double xi = p.x() + q*v.x(),
               yi = p.y() + q*v.y();

      if ( sqr(xi/xSemiAxis) + sqr(yi/ySemiAxis) <= sqr( zheight - zTopCut ) )
      {
        return (sigz > -halfCarTol) ? q : 0;
      }
      else if (xi/(xSemiAxis*xSemiAxis)*v.x()
             + yi/(ySemiAxis*ySemiAxis)*v.y() >= 0)
      {
        //        return kInfinity;
      }
    }
  }


#if 0

  // check to see if Z plane is relevant
  //
  if (p.z() < -zTopCut - halfCarTol)
  {
    if (v.z() <= 0.0)
      return distMin; 

    G4double lambda = (-zTopCut - p.z())/v.z();
    
    if ( sqr((lambda*v.x()+p.x())/xSemiAxis) + 
         sqr((lambda*v.y()+p.y())/ySemiAxis) <=
         sqr(zTopCut + zheight + halfCarTol) ) 
    { 
      return distMin = std::fabs(lambda);    
    }
  }

  if (p.z() > zTopCut + halfCarTol) 
  {
    if (v.z() >= 0.0)
      { return distMin; }

    G4double lambda  = (zTopCut - p.z()) / v.z();

    if ( sqr((lambda*v.x() + p.x())/xSemiAxis) + 
         sqr((lambda*v.y() + p.y())/ySemiAxis) <=
         sqr(zheight - zTopCut + halfCarTol) )
      {
        return distMin = std::fabs(lambda);
      }
  }
  
  if (p.z() > zTopCut - halfCarTol
   && p.z() < zTopCut + halfCarTol )
  {
    if (v.z() > 0.) 
      { return kInfinity; }

    return distMin = 0.;
  }
  
  if (p.z() < -zTopCut + halfCarTol
   && p.z() > -zTopCut - halfCarTol)
  {
    if (v.z() < 0.)
      { return distMin = kInfinity; }
    
    return distMin = 0.;
  }
  
#endif

  // if we are here then it either intersects or grazes the curved surface 
  // or it does not intersect at all
  //
  G4double A = sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) - sqr(v.z());
  G4double B = 2*(v.x()*p.x()/sqr(xSemiAxis) + 
                  v.y()*p.y()/sqr(ySemiAxis) + v.z()*(zheight-p.z()));
  G4double C = sqr(p.x()/xSemiAxis) + sqr(p.y()/ySemiAxis) - 
               sqr(zheight - p.z());
 
  G4double discr = B*B - 4.*A*C;
   
  // if the discriminant is negative it never hits the curved object
  //
  if ( discr < -halfCarTol )
    { return distMin; }
  
  // case below is when it hits or grazes the surface
  //
  if ( (discr >= -halfCarTol ) && (discr < halfCarTol ) )
  {
    return distMin = std::fabs(-B/(2.*A)); 
  }
  
  G4double plus  = (-B+std::sqrt(discr))/(2.*A);
  G4double minus = (-B-std::sqrt(discr))/(2.*A);
 
  // Special case::Point on Surface, Check norm.dot(v)

  if ( ( std::fabs(plus) < halfCarTol )||( std::fabs(minus) < halfCarTol ) )
  {
    G4ThreeVector truenorm(p.x()/(xSemiAxis*xSemiAxis),
                           p.y()/(ySemiAxis*ySemiAxis),
                           -( p.z() - zheight ));
    if ( truenorm*v >= 0)  //  going outside the solid from surface
    {
      return kInfinity;
    }
    else
    {
      return 0;
    }
  }

  // G4double lambda = std::fabs(plus) < std::fabs(minus) ? plus : minus;  
  G4double lambda = 0;

  if ( minus > halfCarTol && minus < distMin ) 
  {
    lambda = minus ;
    // check normal vector   n * v < 0
    G4ThreeVector pin = p + lambda*v;
    if(std::fabs(pin.z())< zTopCut + halfCarTol)
    {
      G4ThreeVector truenorm(pin.x()/(xSemiAxis*xSemiAxis),
                             pin.y()/(ySemiAxis*ySemiAxis),
                             - ( pin.z() - zheight ));
      if ( truenorm*v < 0)
      {   // yes, going inside the solid
        distMin = lambda;
      }
    }
  }
  if ( plus > halfCarTol  && plus < distMin )
  {
    lambda = plus ;
    // check normal vector   n * v < 0
    G4ThreeVector pin = p + lambda*v;
    if(std::fabs(pin.z()) < zTopCut + halfCarTol)
    {
      G4ThreeVector truenorm(pin.x()/(xSemiAxis*xSemiAxis),
                             pin.y()/(ySemiAxis*ySemiAxis),
                             - ( pin.z() - zheight ) );
      if ( truenorm*v < 0)
      {   // yes, going inside the solid
        distMin = lambda;
      }
    }
  }
  if (distMin < halfCarTol) distMin=0.;
  return distMin ;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// Return 0 if point inside

G4double G4EllipticalCone::DistanceToIn(const G4ThreeVector& p) const
{
  G4double hp = std::sqrt(p.x()*p.x()*invXX + p.y()*p.y()*invYY) + p.z();
  G4double ds = (hp - zheight)*cosAxisMin;
  G4double dz = std::abs(p.z()) - zTopCut;
  G4double dist = std::max(ds,dz);
  return (dist > 0) ? dist : 0.;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside',
// allowing for tolerance

G4double G4EllipticalCone::DistanceToOut(const G4ThreeVector& p,
                                         const G4ThreeVector& v,
                                         const G4bool calcNorm,
                                               G4bool* validNorm,
                                               G4ThreeVector* n  ) const
{
  G4double distMin, lambda;
  enum surface_e {kPlaneSurf, kCurvedSurf, kNoSurf} surface;
  
  distMin = kInfinity;
  surface = kNoSurf;

  if (v.z() < 0.0)
  {
    lambda = (-p.z() - zTopCut)/v.z();

    if ( (sqr((p.x() + lambda*v.x())/xSemiAxis) + 
          sqr((p.y() + lambda*v.y())/ySemiAxis)) < 
          sqr(zheight + zTopCut + halfCarTol) )
    {
      distMin = std::fabs(lambda);

      if (!calcNorm) { return distMin; }
    } 
    distMin = std::fabs(lambda);
    surface = kPlaneSurf;
  }

  if (v.z() > 0.0)
  {
    lambda = (zTopCut - p.z()) / v.z();

    if ( (sqr((p.x() + lambda*v.x())/xSemiAxis)
        + sqr((p.y() + lambda*v.y())/ySemiAxis) )
       < (sqr(zheight - zTopCut + halfCarTol)) )
    {
      distMin = std::fabs(lambda);
      if (!calcNorm) { return distMin; }
    }
    distMin = std::fabs(lambda);
    surface = kPlaneSurf;
  }
  
  // if we are here then it either intersects or grazes the 
  // curved surface...
  //
  G4double A = sqr(v.x()/xSemiAxis) + sqr(v.y()/ySemiAxis) - sqr(v.z());
  G4double B = 2.*(v.x()*p.x()/sqr(xSemiAxis) +  
                   v.y()*p.y()/sqr(ySemiAxis) + v.z()*(zheight-p.z()));
  G4double C = sqr(p.x()/xSemiAxis) + sqr(p.y()/ySemiAxis)
             - sqr(zheight - p.z());
 
  G4double discr = B*B - 4.*A*C;
  
  if ( discr >= - halfCarTol && discr < halfCarTol )
  { 
    if(!calcNorm) { return distMin = std::fabs(-B/(2.*A)); }
  }

  else if ( discr > halfCarTol )
  {
    G4double plus  = (-B+std::sqrt(discr))/(2.*A);
    G4double minus = (-B-std::sqrt(discr))/(2.*A);

    if ( plus > halfCarTol && minus > halfCarTol )
    {
      // take the shorter distance
      //
      lambda   = std::fabs(plus) < std::fabs(minus) ? plus : minus;
    }
    else
    {
      // at least one solution is close to zero or negative
      // so, take small positive solution or zero 
      //
      lambda   = plus > -halfCarTol ? plus : 0;
    }

    if ( std::fabs(lambda) < distMin )
    {
      if( std::fabs(lambda) > halfCarTol)
      {
        distMin  = std::fabs(lambda);
        surface  = kCurvedSurf;
      }
      else  // Point is On the Surface, Check Normal
      {
        G4ThreeVector truenorm(p.x()/(xSemiAxis*xSemiAxis),
                               p.y()/(ySemiAxis*ySemiAxis),
                               -( p.z() - zheight ));
        if( truenorm.dot(v) > 0 )
        {
          distMin  = 0.0;
          surface  = kCurvedSurf;
        }
      } 
    }
  }

  // set normal if requested
  //
  if (calcNorm)
  {
    if (surface == kNoSurf)
    {
      *validNorm = false;
    }
    else
    {
      *validNorm = true;
      switch (surface)
      {
        case kPlaneSurf:
        {
          *n = G4ThreeVector(0.,0.,(v.z() > 0.0 ? 1. : -1.));
        }
        break;

        case kCurvedSurf:
        {
          G4ThreeVector pexit = p + distMin*v;
          G4ThreeVector truenorm( pexit.x()/(xSemiAxis*xSemiAxis),
                                  pexit.y()/(ySemiAxis*ySemiAxis),
                                  -( pexit.z() - zheight ) );
          truenorm /= truenorm.mag();
          *n= truenorm;
        } 
        break;

        default:            // Should never reach this case ...
          DumpInfo();
          std::ostringstream message;
          G4long oldprc = message.precision(16);
          message << "Undefined side for valid surface normal to solid."
                  << G4endl
                  << "Position:"  << G4endl
                  << "   p.x() = "   << p.x()/mm << " mm" << G4endl
                  << "   p.y() = "   << p.y()/mm << " mm" << G4endl
                  << "   p.z() = "   << p.z()/mm << " mm" << G4endl
                  << "Direction:" << G4endl
                  << "   v.x() = "   << v.x() << G4endl
                  << "   v.y() = "   << v.y() << G4endl
                  << "   v.z() = "   << v.z() << G4endl
                  << "Proposed distance :" << G4endl
                  << "   distMin = "    << distMin/mm << " mm";
          message.precision(oldprc);
          G4Exception("G4EllipticalCone::DistanceToOut(p,v,..)",
                      "GeomSolids1002", JustWarning, message);
          break;
      }
    }
  }

  if (distMin < halfCarTol) { distMin=0; }

  return distMin;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4EllipticalCone::DistanceToOut(const G4ThreeVector& p) const
{
#ifdef G4SPECSDEBUG
  if( Inside(p) == kOutside )
  {
     std::ostringstream message;
     G4long oldprc = message.precision(16);
     message << "Point p is outside (!?) of solid: " << GetName() << "\n"
             << "Position:\n"
             << "   p.x() = "  << p.x()/mm << " mm\n"
             << "   p.y() = "  << p.y()/mm << " mm\n"
             << "   p.z() = "  << p.z()/mm << " mm";
     message.precision(oldprc) ;
     G4Exception("G4Ellipsoid::DistanceToOut(p)", "GeomSolids1002",
                 JustWarning, message);
     DumpInfo();
  }
#endif
  G4double hp = std::sqrt(p.x()*p.x()*invXX + p.y()*p.y()*invYY) + p.z();
  G4double ds = (zheight - hp)*cosAxisMin;
  G4double dz = zTopCut - std::abs(p.z());
  G4double dist = std::min(ds,dz);
  return (dist > 0) ? dist : 0.;
}

/////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4EllipticalCone::GetEntityType() const
{
  return G4String("G4EllipticalCone");
}

/////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4EllipticalCone::Clone() const
{
  return new G4EllipticalCone(*this);
}

/////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4EllipticalCone::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4EllipticalCone\n"
     << " Parameters: \n"

     << "    semi-axis x: " << xSemiAxis/mm << " mm \n"
     << "    semi-axis y: " << ySemiAxis/mm << " mm \n"
     << "    height    z: " << zheight/mm << " mm \n"
     << "    half length in  z: " << zTopCut/mm << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// Return random point on the surface of the solid

G4ThreeVector G4EllipticalCone::GetPointOnSurface() const
{
  G4double x0 = xSemiAxis*zheight; // x semi axis at z=0
  G4double y0 = ySemiAxis*zheight; // y semi axis at z=0
  G4double s0 = G4GeomTools::EllipticConeLateralArea(x0,y0,zheight);
  G4double kmin = (zTopCut >= zheight ) ? 0. : (zheight - zTopCut)/zheight;
  G4double kmax = (zTopCut >= zheight ) ? 2. : (zheight + zTopCut)/zheight;

  // Set areas (base at -Z, side surface, base at +Z)
  //
  G4double szmin =  pi*x0*y0*kmax*kmax;
  G4double szmax =  pi*x0*y0*kmin*kmin;
  G4double sside =  s0*(kmax*kmax - kmin*kmin);
  G4double ssurf[3] = { szmin, sside, szmax };
  for (auto i=1; i<3; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[2]*G4UniformRand();
  G4int k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Pick random point on selected surface
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // base at -Z, uniform distribution, rejection sampling
    {
      G4double zh = zheight + zTopCut;
      G4TwoVector rho = G4RandomPointInEllipse(zh*xSemiAxis,zh*ySemiAxis);
      p.set(rho.x(),rho.y(),-zTopCut);
      break;
    }
    case 1: // side surface, uniform distribution, rejection sampling
    {
      G4double zh = G4RandomRadiusInRing(zheight-zTopCut, zheight+zTopCut);
      G4double a = x0;
      G4double b = y0;

      G4double hh = zheight*zheight;
      G4double aa = a*a;
      G4double bb = b*b;
      G4double R  = std::max(a,b);
      G4double mu_max = R*std::sqrt(hh + R*R);

      G4double x,y;
      for (auto i=0; i<1000; ++i)
      {
	G4double phi = CLHEP::twopi*G4UniformRand();
        x = std::cos(phi);
        y = std::sin(phi);
        G4double xx = x*x;
        G4double yy = y*y;
        G4double E = hh + aa*xx + bb*yy;
        G4double F = (aa-bb)*x*y;
        G4double G = aa*yy + bb*xx;
        G4double mu = std::sqrt(E*G - F*F);
        if (mu_max*G4UniformRand() <= mu) break;
      }
      p.set(zh*xSemiAxis*x,zh*ySemiAxis*y,zheight-zh);
      break;
    }
    case 2: // base at +Z, uniform distribution, rejection sampling
    {
      G4double zh = zheight - zTopCut;
      G4TwoVector rho = G4RandomPointInEllipse(zh*xSemiAxis,zh*ySemiAxis);
      p.set(rho.x(),rho.y(),zTopCut);
      break;
    }
  }
  return p;
}

/////////////////////////////////////////////////////////////////////////
//
// Get cubic volume

G4double G4EllipticalCone::GetCubicVolume()
{
  if (fCubicVolume == 0.0)
  {
    G4double x0 = xSemiAxis*zheight; // x semi axis at z=0
    G4double y0 = ySemiAxis*zheight; // y semi axis at z=0
    G4double v0 = CLHEP::pi*x0*y0*zheight/3.;
    G4double kmin = (zTopCut >= zheight ) ? 0. : (zheight - zTopCut)/zheight;
    G4double kmax = (zTopCut >= zheight ) ? 2. : (zheight + zTopCut)/zheight;
    fCubicVolume = (kmax - kmin)*(kmax*kmax + kmax*kmin + kmin*kmin)*v0;
  }
  return fCubicVolume;
}

/////////////////////////////////////////////////////////////////////////
//
// Get surface area

G4double G4EllipticalCone::GetSurfaceArea()
{
  if (fSurfaceArea == 0.0)
  {
    G4double x0 = xSemiAxis*zheight; // x semi axis at z=0
    G4double y0 = ySemiAxis*zheight; // y semi axis at z=0
    G4double s0 = G4GeomTools::EllipticConeLateralArea(x0,y0,zheight);
    G4double kmin = (zTopCut >= zheight ) ? 0. : (zheight - zTopCut)/zheight;
    G4double kmax = (zTopCut >= zheight ) ? 2. : (zheight + zTopCut)/zheight;
    fSurfaceArea = (kmax - kmin)*(kmax + kmin)*s0
                 + CLHEP::pi*x0*y0*(kmin*kmin + kmax*kmax);
  }
  return fSurfaceArea;
}

/////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4EllipticalCone::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4VisExtent G4EllipticalCone::GetExtent() const
{
  // Define the sides of the box into which the solid instance would fit.
  //
  G4ThreeVector pmin,pmax;
  BoundingLimits(pmin,pmax);
  return G4VisExtent(pmin.x(),pmax.x(),
                     pmin.y(),pmax.y(),
                     pmin.z(),pmax.z());
}

G4Polyhedron* G4EllipticalCone::CreatePolyhedron () const
{
  return new G4PolyhedronEllipticalCone(xSemiAxis, ySemiAxis, zheight, zTopCut);
}

G4Polyhedron* G4EllipticalCone::GetPolyhedron () const
{
  if ( (fpPolyhedron == nullptr)
    || fRebuildPolyhedron
    || (fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
        fpPolyhedron->GetNumberOfRotationSteps()) )
    {
      G4AutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }
  return fpPolyhedron;
}

#endif // !defined(G4GEOM_USE_UELLIPTICALCONE) || !defined(G4GEOM_USE_SYS_USOLIDS)
