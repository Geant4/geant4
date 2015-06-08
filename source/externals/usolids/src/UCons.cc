//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UCons
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UCons.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull, kRMin, kRMax, kSPhi, kEPhi, kPZ, kMZ};

// used by normal

enum ENorm {kNRMin, kNRMax, kNSPhi, kNEPhi, kNZ};

//////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//               - note if pDPhi>2PI then reset to 2PI

UCons::UCons(const std::string& pName,
             double  pRmin1, double pRmax1,
             double  pRmin2, double pRmax2,
             double pDz,
             double pSPhi, double pDPhi)
  : VUSolid(pName.c_str()), fRmin1(pRmin1), fRmin2(pRmin2),
    fRmax1(pRmax1), fRmax2(pRmax2), fDz(pDz), fSPhi(0.), fDPhi(0.)
{
  kRadTolerance = frTolerance;
  kAngTolerance = faTolerance;
  // Check z-len
  //
  if (pDz < 0)
  {
    std::ostringstream message;
    message << "Invalid Z half-length for Solid: " << GetName() << std::endl
            << "  hZ = " << pDz;
    UUtils::Exception("UCons::UCons()", "GeomSolids0002", UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Check radii
  //
  if (((pRmin1 >= pRmax1) || (pRmin2 >= pRmax2) || (pRmin1 < 0)) && (pRmin2 < 0))
  {
    std::ostringstream message;
    message << "Invalid values of radii for Solid: " << GetName() << std::endl
            << "  pRmin1 = " << pRmin1 << ", pRmin2 = " << pRmin2
            << ", pRmax1 = " << pRmax1 << ", pRmax2 = " << pRmax2;
    UUtils::Exception("UCons::UCons()", "GeomSolids0002", UFatalErrorInArguments, 1, message.str().c_str());

  }
  if ((pRmin1 == 0.0) && (pRmin2 > 0.0))
  {
    fRmin1 = 1e3 * kRadTolerance;
  }
  if ((pRmin2 == 0.0) && (pRmin1 > 0.0))
  {
    fRmin2 = 1e3 * kRadTolerance;
  }

  // Check angles
  //
  CheckPhiAngles(pSPhi, pDPhi);

  Initialize();
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UCons::UCons(/* __void__& a */)
  : VUSolid(""), kRadTolerance(0.), kAngTolerance(0.),
    fRmin1(0.), fRmin2(0.), fRmax1(0.), fRmax2(0.), fDz(0.),
    fSPhi(0.), fDPhi(0.), sinCPhi(0.), cosCPhi(0.), cosHDPhiOT(0.),
    cosHDPhiIT(0.), sinSPhi(0.), cosSPhi(0.), sinEPhi(0.), cosEPhi(0.),
    fPhiFullCone(false)
{
  Initialize();
}

///////////////////////////////////////////////////////////////////////
//
// Destructor

UCons::~UCons()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

UCons::UCons(const UCons& rhs)
  : VUSolid(rhs), kRadTolerance(rhs.kRadTolerance),
    kAngTolerance(rhs.kAngTolerance), fRmin1(rhs.fRmin1), fRmin2(rhs.fRmin2),
    fRmax1(rhs.fRmax1), fRmax2(rhs.fRmax2), fDz(rhs.fDz), fSPhi(rhs.fSPhi),
    fDPhi(rhs.fDPhi), sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi), sinEPhi(rhs.sinEPhi),
    cosEPhi(rhs.cosEPhi), fPhiFullCone(rhs.fPhiFullCone)
{
  Initialize();
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

UCons& UCons::operator = (const UCons& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  VUSolid::operator=(rhs);

  // Copy data
  //
  kRadTolerance = rhs.kRadTolerance;
  kAngTolerance = rhs.kAngTolerance;
  fRmin1 = rhs.fRmin1;
  fRmin2 = rhs.fRmin2;
  fRmax1 = rhs.fRmax1;
  fRmax2 = rhs.fRmax2;
  fDz = rhs.fDz;
  fSPhi = rhs.fSPhi;
  fDPhi = rhs.fDPhi;
  sinCPhi = rhs.sinCPhi;
  cosCPhi = rhs.cosCPhi;
  cosHDPhiOT = rhs.cosHDPhiOT;
  cosHDPhiIT = rhs.cosHDPhiIT;
  sinSPhi = rhs.sinSPhi;
  cosSPhi = rhs.cosSPhi;
  sinEPhi = rhs.sinEPhi;
  cosEPhi = rhs.cosEPhi;
  fPhiFullCone = rhs.fPhiFullCone;

  Initialize();
  return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Return Unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

bool UCons::Normal(const UVector3& p, UVector3& n) const
{
  int noSurfaces = 0;
  double rho, pPhi;
  double distZ, distRMin, distRMax;
  double distSPhi = UUtils::kInfinity, distEPhi = UUtils::kInfinity;
  double pRMin, widRMin;
  double pRMax, widRMax;

  static const double delta = 0.5 * VUSolid::Tolerance();
  static const double dAngle = 0.5 * kAngTolerance;

  UVector3 norm, sumnorm(0., 0., 0.), nZ = UVector3(0., 0., 1.);
  UVector3 nR, nr(0., 0., 0.), nPs, nPe;

  distZ = std::fabs(std::fabs(p.z()) - fDz);
  rho  = std::sqrt(p.x() * p.x() + p.y() * p.y());

  pRMin   = rho - p.z() * tanRMin;
  widRMin = fRmin2 - fDz * tanRMin;
  distRMin = std::fabs(pRMin - widRMin) / secRMin;

  pRMax   = rho - p.z() * tanRMax;
  widRMax = fRmax2 - fDz * tanRMax;
  distRMax = std::fabs(pRMax - widRMax) / secRMax;

  if (!fPhiFullCone)   // Protected against (0,0,z)
  {
    if (rho)
    {
      pPhi = std::atan2(p.y(), p.x());

      if (pPhi  < fSPhi - delta)
      {
        pPhi += 2 * UUtils::kPi;
      }
      else if (pPhi > fSPhi + fDPhi + delta)
      {
        pPhi -= 2 * UUtils::kPi;
      }

      distSPhi = std::fabs(pPhi - fSPhi);
      distEPhi = std::fabs(pPhi - fSPhi - fDPhi);
    }
    else if (!(fRmin1) || !(fRmin2))
    {
      distSPhi = 0.;
      distEPhi = 0.;
    }
    nPs = UVector3(std::sin(fSPhi), -std::cos(fSPhi), 0);
    nPe = UVector3(-std::sin(fSPhi + fDPhi), std::cos(fSPhi + fDPhi), 0);
  }
  if (rho > delta)
  {
    nR = UVector3(p.x() / rho / secRMax, p.y() / rho / secRMax, -tanRMax / secRMax);
    if (fRmin1 || fRmin2)
    {
      nr = UVector3(-p.x() / rho / secRMin, -p.y() / rho / secRMin, tanRMin / secRMin);
    }
  }

  if (distRMax <= delta)
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if ((fRmin1 || fRmin2) && (distRMin <= delta))
  {
    noSurfaces ++;
    sumnorm += nr;
  }
  if (!fPhiFullCone)
  {
    if (distSPhi <= dAngle)
    {
      noSurfaces ++;
      sumnorm += nPs;
    }
    if (distEPhi <= dAngle)
    {
      noSurfaces ++;
      sumnorm += nPe;
    }
  }
  if (distZ <= delta)
  {
    noSurfaces ++;
    if (p.z() >= 0.)
    {
      sumnorm += nZ;
    }
    else
    {
      sumnorm -= nZ;
    }
  }
  if (noSurfaces == 0)
  {
#ifdef UDEBUG

    UUtils::Exception("UCons::SurfaceNormal(p)", "GeomSolids1002",
                      UWarning, 1, "Point p is not on surface !?");
#endif
    norm = ApproxSurfaceNormal(p);
  }
  else if (noSurfaces == 1)
  {
    norm = sumnorm;
  }
  else
  {
    norm = sumnorm.Unit();
  }

  n = norm;

  return (bool) noSurfaces;
}

////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

UVector3 UCons::ApproxSurfaceNormal(const UVector3& p) const
{
  ENorm side;
  UVector3 norm;
  double rho, phi;
  double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin;
  double pRMin, widRMin;
  double pRMax, widRMax;

  distZ = std::fabs(std::fabs(p.z()) - fDz);
  rho  = std::sqrt(p.x() * p.x() + p.y() * p.y());

  pRMin   = rho - p.z() * tanRMin;
  widRMin = fRmin2 - fDz * tanRMin;
  distRMin = std::fabs(pRMin - widRMin) / secRMin;

  pRMax   = rho - p.z() * tanRMax;
  widRMax = fRmax2 - fDz * tanRMax;
  distRMax = std::fabs(pRMax - widRMax) / secRMax;

  if (distRMin < distRMax)  // First minimum
  {
    if (distZ < distRMin)
    {
      distMin = distZ;
      side    = kNZ;
    }
    else
    {
      distMin = distRMin;
      side    = kNRMin;
    }
  }
  else
  {
    if (distZ < distRMax)
    {
      distMin = distZ;
      side    = kNZ;
    }
    else
    {
      distMin = distRMax;
      side    = kNRMax;
    }
  }
  if (!fPhiFullCone && rho)    // Protected against (0,0,z)
  {
    phi = std::atan2(p.y(), p.x());

    if (phi < 0)
    {
      phi += 2 * UUtils::kPi;
    }

    if (fSPhi < 0)
    {
      distSPhi = std::fabs(phi - (fSPhi + 2 * UUtils::kPi)) * rho;
    }
    else
    {
      distSPhi = std::fabs(phi - fSPhi) * rho;
    }

    distEPhi = std::fabs(phi - fSPhi - fDPhi) * rho;

    // Find new minimum

    if (distSPhi < distEPhi)
    {
      if (distSPhi < distMin)
      {
        side = kNSPhi;
      }
    }
    else
    {
      if (distEPhi < distMin)
      {
        side = kNEPhi;
      }
    }
  }
  switch (side)
  {
    case kNRMin:      // Inner radius
      rho *= secRMin;
      norm = UVector3(-p.x() / rho, -p.y() / rho, tanRMin / secRMin);
      break;
    case kNRMax:      // Outer radius
      rho *= secRMax;
      norm = UVector3(p.x() / rho, p.y() / rho, -tanRMax / secRMax);
      break;
    case kNZ:        // +/- dz
      if (p.z() > 0)
      {
        norm = UVector3(0, 0, 1);
      }
      else
      {
        norm = UVector3(0, 0, -1);
      }
      break;
    case kNSPhi:
      norm = UVector3(std::sin(fSPhi), -std::cos(fSPhi), 0);
      break;
    case kNEPhi:
      norm = UVector3(-std::sin(fSPhi + fDPhi), std::cos(fSPhi + fDPhi), 0);
      break;
    default:          // Should never reach this case...
      UUtils::Exception("UCons::ApproxSurfaceNormal()",
                        "GeomSolids1002", UWarning, 1,
                        "Undefined side for valid surface normal to solid.");
      break;
  }
  return norm;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return UUtils::kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes
//        - if at valid r, phi, return
//
// -> If point is outside cone, compute intersection with rmax1*0.5
//        - if at valid phi,z return
//        - if inside outer cone, handle case when on tolerant outer cone
//          boundary and heading inwards(->0 to in)
//
// -> Compute intersection with inner cone, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - `if valid' implies tolerant checking of intersection points
// - z, phi intersection from Tubs

double UCons::DistanceToIn(const UVector3& p,
                           const UVector3& v, double /* aPstep */) const
{
  double snxt = UUtils::kInfinity;
  const double dRmax = 100 * std::max(fRmax1, fRmax2);
  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfRadTolerance = kRadTolerance * 0.5;

  double rMaxAv, rMaxOAv; // Data for cones
  double rMinAv, rMinOAv;
  double rout, rin;

  double tolORMin, tolORMin2, tolIRMin, tolIRMin2; // `generous' radii squared
  double tolORMax2, tolIRMax, tolIRMax2;
  double tolODz, tolIDz;

  double Dist, sd, xi, yi, zi, ri = 0., risec, rhoi2, cosPsi; // Intersection point vars

  double t1, t2, t3, b, c, d; // Quadratic solver variables
  double nt1, nt2, nt3;
  double Comp;

  UVector3 norm;

  // Cone Precalcs
  rMinAv  = (fRmin1 + fRmin2) * 0.5;

  if (rMinAv > halfRadTolerance)
  {
    rMinOAv = rMinAv - halfRadTolerance;
  }
  else
  {
    rMinOAv = 0.0;
  }
  rMaxAv  = (fRmax1 + fRmax2) * 0.5;
  rMaxOAv = rMaxAv + halfRadTolerance;

  // Intersection with z-surfaces

  tolIDz = fDz - halfCarTolerance;
  tolODz = fDz + halfCarTolerance;

  if (std::fabs(p.z()) >= tolIDz)
  {
    if (p.z() * v.z() < 0)   // at +Z going in -Z or visa versa
    {
      sd = (std::fabs(p.z()) - fDz) / std::fabs(v.z()); // Z intersect distance

      if (sd < 0.0)
      {
        sd = 0.0;  // negative dist -> zero
      }

      xi   = p.x() + sd * v.x(); // Intersection coords
      yi   = p.y() + sd * v.y();
      rhoi2 = xi * xi + yi * yi ;

      // Check validity of intersection
      // Calculate (outer) tolerant radi^2 at intersecion

      if (v.z() > 0)
      {
        tolORMin  = fRmin1 - halfRadTolerance * secRMin;
        tolIRMin  = fRmin1 + halfRadTolerance * secRMin;
        tolIRMax  = fRmax1 - halfRadTolerance * secRMin;
        // tolORMax2 = (fRmax1 + halfRadTolerance * secRMax) *
        //             (fRmax1 + halfRadTolerance * secRMax);
      }
      else
      {
        tolORMin  = fRmin2 - halfRadTolerance * secRMin;
        tolIRMin  = fRmin2 + halfRadTolerance * secRMin;
        tolIRMax  = fRmax2 - halfRadTolerance * secRMin;
        // tolORMax2 = (fRmax2 + halfRadTolerance * secRMax) *
        //             (fRmax2 + halfRadTolerance * secRMax);
      }
      if (tolORMin > 0)
      {
        // tolORMin2 = tolORMin * tolORMin;
        tolIRMin2 = tolIRMin * tolIRMin;
      }
      else
      {
        // tolORMin2 = 0.0;
        tolIRMin2 = 0.0;
      }
      if (tolIRMax > 0)
      {
        tolIRMax2 = tolIRMax * tolIRMax;
      }
      else
      {
        tolIRMax2 = 0.0;
      }

      if ((tolIRMin2 <= rhoi2) && (rhoi2 <= tolIRMax2))
      {
        if (!fPhiFullCone && rhoi2)
        {
          // Psi = angle made with central (average) phi of shape

          cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);

          if (cosPsi >= cosHDPhiIT)
          {
            return sd;
          }
        }
        else
        {
          return sd;
        }
      }
    }
    else  // On/outside extent, and heading away  -> cannot intersect
    {
      return snxt;
    }
  }

// ----> Can not intersect z surfaces


// Intersection with outer cone (possible return) and
//                   inner cone (must also check phi)
//
// Intersection point (xi,yi,zi) on line x=p.x()+t*v.x() etc.
//
// Intersects with x^2+y^2=(a*z+b)^2
//
// where a=tanRMax or tanRMin
//       b=rMaxAv or rMinAv
//
// (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0;
//     t1                       t2                      t3
//
//  \--------u-------/       \-----------v----------/ \---------w--------/
//

  t1   = 1.0 - v.z() * v.z();
  t2   = p.x() * v.x() + p.y() * v.y();
  t3   = p.x() * p.x() + p.y() * p.y();
  rin = tanRMin * p.z() + rMinAv;
  rout = tanRMax * p.z() + rMaxAv;

  // Outer Cone Intersection
  // Must be outside/on outer cone for valid intersection

  nt1 = t1 - (tanRMax * v.z()) * (tanRMax * v.z());
  nt2 = t2 - tanRMax * v.z() * rout;
  nt3 = t3 - rout * rout;

  if (std::fabs(nt1) > kRadTolerance) // Equation quadratic => 2 roots
  {
    b = nt2 / nt1;
    c = nt3 / nt1;
    d = b * b - c ;
    if ((nt3 > rout * rout * kRadTolerance * kRadTolerance * secRMax * secRMax)
        || (rout < 0))
    {
      // If outside real cone (should be rho-rout>kRadTolerance*0.5
      // NOT rho^2 etc) saves a std::sqrt() at expense of accuracy

      if (d >= 0)
      {

        if ((rout < 0) && (nt3 <= 0))
        {
          // Inside `shadow cone' with -ve radius
          // -> 2nd root could be on real cone

          if (b > 0)
          {
            sd = c / (-b - std::sqrt(d));
          }
          else
          {
            sd = -b + std::sqrt(d);
          }
        }
        else
        {
          if ((b <= 0) && (c >= 0)) // both >=0, try smaller root
          {
            sd = c / (-b + std::sqrt(d));
          }
          else
          {
            if (c <= 0)   // second >=0
            {
              sd = -b + std::sqrt(d);
            }
            else  // both negative, travel away
            {
              return UUtils::kInfinity;
            }
          }
        }
        if (sd >= 0)    // If 'forwards'. Check z intersection
        {
          if (sd > dRmax) // Avoid rounding errors due to precision issues on
          {
            // 64 bits systems. Split long distances and recompute
            double fTerm = sd - std::fmod(sd, dRmax);
            sd = fTerm + DistanceToIn(p + fTerm * v, v);
          }
          zi = p.z() + sd * v.z();

          if (std::fabs(zi) <= tolODz)
          {
            // Z ok. Check phi intersection if reqd

            if (fPhiFullCone)
            {
              return sd;
            }
            else
            {
              xi     = p.x() + sd * v.x();
              yi     = p.y() + sd * v.y();
              ri     = rMaxAv + zi * tanRMax;
              cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

              if (cosPsi >= cosHDPhiIT)
              {
                return sd;
              }
            }
          }
        }               // end if (sd>0)
      }
    }
    else
    {
      // Inside outer cone
      // check not inside, and heading through UCons (-> 0 to in)

      if ((t3 > (rin + halfRadTolerance * secRMin)*
           (rin + halfRadTolerance * secRMin))
          && (nt2 < 0) && (d >= 0) && (std::fabs(p.z()) <= tolIDz))
      {
        // Inside cones, delta r -ve, inside z extent
        // Point is on the Surface => check Direction using Normal.Dot(v)

        xi     = p.x();
        yi     = p.y()  ;
        risec = std::sqrt(xi * xi + yi * yi) * secRMax;
        norm = UVector3(xi / risec, yi / risec, -tanRMax / secRMax);
        if (!fPhiFullCone)
        {
          cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(t3);
          if (cosPsi >= cosHDPhiIT)
          {
            if (norm.Dot(v) <= 0)
            {
              return 0.0;
            }
          }
        }
        else
        {
          if (norm.Dot(v) <= 0)
          {
            return 0.0;
          }
        }
      }
    }
  }
  else  //  Single root case
  {
    if (std::fabs(nt2) > kRadTolerance)
    {
      sd = -0.5 * nt3 / nt2;

      if (sd < 0)
      {
        return UUtils::kInfinity;  // travel away
      }
      else  // sd >= 0, If 'forwards'. Check z intersection
      {
        zi = p.z() + sd * v.z();

        if ((std::fabs(zi) <= tolODz) && (nt2 < 0))
        {
          // Z ok. Check phi intersection if reqd

          if (fPhiFullCone)
          {
            return sd;
          }
          else
          {
            xi     = p.x() + sd * v.x();
            yi     = p.y() + sd * v.y();
            ri     = rMaxAv + zi * tanRMax;
            cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

            if (cosPsi >= cosHDPhiIT)
            {
              return sd;
            }
          }
        }
      }
    }
    else  //    travel || cone surface from its origin
    {
      sd = UUtils::kInfinity;
    }
  }

  // Inner Cone Intersection
  // o Space is divided into 3 areas:
  //   1) Radius greater than real inner cone & imaginary cone & outside
  //      tolerance
  //   2) Radius less than inner or imaginary cone & outside tolarance
  //   3) Within tolerance of real or imaginary cones
  //      - Extra checks needed for 3's intersections
  //        => lots of duplicated code

  if (rMinAv)
  {
    nt1 = t1 - (tanRMin * v.z()) * (tanRMin * v.z());
    nt2 = t2 - tanRMin * v.z() * rin;
    nt3 = t3 - rin * rin;

    if (nt1)
    {
      if (nt3 > rin * kRadTolerance * secRMin)
      {
        // At radius greater than real & imaginary cones
        // -> 2nd root, with zi check

        b = nt2 / nt1;
        c = nt3 / nt1;
        d = b * b - c;
        if (d >= 0)  // > 0
        {
          if (b > 0)
          {
            sd = c / (-b - std::sqrt(d));
          }
          else
          {
            sd = -b + std::sqrt(d);
          }

          if (sd >= 0)    // > 0
          {
            if (sd > dRmax) // Avoid rounding errors due to precision issues on
            {
              // 64 bits systems. Split long distance and recompute
              double fTerm = sd - std::fmod(sd, dRmax);
              sd = fTerm + DistanceToIn(p + fTerm * v, v);
            }
            zi = p.z() + sd * v.z();

            if (std::fabs(zi) <= tolODz)
            {
              if (!fPhiFullCone)
              {
                xi     = p.x() + sd * v.x();
                yi     = p.y() + sd * v.y();
                ri     = rMinAv + zi * tanRMin;
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

                if (cosPsi >= cosHDPhiIT)
                {
                  if (sd > halfRadTolerance)
                  {
                    snxt = sd;
                  }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                    norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                    if (norm.Dot(v) <= 0)
                    {
                      snxt = sd;
                    }
                  }
                }
              }
              else
              {
                if (sd > halfRadTolerance)
                {
                  return sd;
                }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd * v.x();
                  yi     = p.y() + sd * v.y();
                  risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                  norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                  if (norm.Dot(v) <= 0)
                  {
                    return sd;
                  }
                }
              }
            }
          }
        }
      }
      else  if (nt3 < -rin * kRadTolerance * secRMin)
      {
        // Within radius of inner cone (real or imaginary)
        // -> Try 2nd root, with checking intersection is with real cone
        // -> If check fails, try 1st root, also checking intersection is
        //    on real cone

        b = nt2 / nt1;
        c = nt3 / nt1;
        d = b * b - c;

        if (d >= 0)    // > 0
        {
          if (b > 0)
          {
            sd = c / (-b - std::sqrt(d));
          }
          else
          {
            sd = -b + std::sqrt(d);
          }
          zi = p.z() + sd * v.z();
          ri = rMinAv + zi * tanRMin;

          if (ri > 0)
          {
            if ((sd >= 0) && (std::fabs(zi) <= tolODz))    // sd > 0
            {
              if (sd > dRmax) // Avoid rounding errors due to precision issues
              {
                // seen on 64 bits systems. Split and recompute
                double fTerm = sd - std::fmod(sd, dRmax);
                sd = fTerm + DistanceToIn(p + fTerm * v, v);
              }
              if (!fPhiFullCone)
              {
                xi     = p.x() + sd * v.x();
                yi     = p.y() + sd * v.y();
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

                if (cosPsi >= cosHDPhiOT)
                {
                  if (sd > halfRadTolerance)
                  {
                    snxt = sd;
                  }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                    norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                    if (norm.Dot(v) <= 0)
                    {
                      snxt = sd;
                    }
                  }
                }
              }
              else
              {
                if (sd > halfRadTolerance)
                {
                  return sd;
                }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd * v.x();
                  yi     = p.y() + sd * v.y();
                  risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                  norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                  if (norm.Dot(v) <= 0)
                  {
                    return sd;
                  }
                }
              }
            }
          }
          else
          {
            if (b > 0)
            {
              sd = -b - std::sqrt(d);
            }
            else
            {
              sd = c / (-b + std::sqrt(d));
            }
            zi = p.z() + sd * v.z();
            ri = rMinAv + zi * tanRMin;

            if ((sd >= 0) && (ri > 0) && (std::fabs(zi) <= tolODz))   // sd>0
            {
              if (sd > dRmax) // Avoid rounding errors due to precision issues
              {
                // seen on 64 bits systems. Split and recompute
                double fTerm = sd - std::fmod(sd, dRmax);
                sd = fTerm + DistanceToIn(p + fTerm * v, v);
              }
              if (!fPhiFullCone)
              {
                xi     = p.x() + sd * v.x();
                yi     = p.y() + sd * v.y();
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

                if (cosPsi >= cosHDPhiIT)
                {
                  if (sd > halfRadTolerance)
                  {
                    snxt = sd;
                  }
                  else
                  {
                    // Calculate a normal vector in order to check Direction

                    risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                    norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                    if (norm.Dot(v) <= 0)
                    {
                      snxt = sd;
                    }
                  }
                }
              }
              else
              {
                if (sd > halfRadTolerance)
                {
                  return sd;
                }
                else
                {
                  // Calculate a normal vector in order to check Direction

                  xi     = p.x() + sd * v.x();
                  yi     = p.y() + sd * v.y();
                  risec = std::sqrt(xi * xi + yi * yi) * secRMin;
                  norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
                  if (norm.Dot(v) <= 0)
                  {
                    return sd;
                  }
                }
              }
            }
          }
        }
      }
      else
      {
        // Within kRadTol*0.5 of inner cone (real OR imaginary)
        // ----> Check not travelling through (=>0 to in)
        // ----> if not:
        //    -2nd root with validity check

        if (std::fabs(p.z()) <= tolODz)
        {
          if (nt2 > 0)
          {
            // Inside inner real cone, heading outwards, inside z range

            if (!fPhiFullCone)
            {
              cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(t3);

              if (cosPsi >= cosHDPhiIT)
              {
                return 0.0;
              }
            }
            else
            {
              return 0.0;
            }
          }
          else
          {
            // Within z extent, but not travelling through
            // -> 2nd root or UUtils::kInfinity if 1st root on imaginary cone

            b = nt2 / nt1;
            c = nt3 / nt1;
            d = b * b - c;

            if (d >= 0)     // > 0
            {
              if (b > 0)
              {
                sd = -b - std::sqrt(d);
              }
              else
              {
                sd = c / (-b + std::sqrt(d));
              }
              zi = p.z() + sd * v.z();
              ri = rMinAv + zi * tanRMin;

              if (ri > 0)     // 2nd root
              {
                if (b > 0)
                {
                  sd = c / (-b - std::sqrt(d));
                }
                else
                {
                  sd = -b + std::sqrt(d);
                }

                zi = p.z() + sd * v.z();

                if ((sd >= 0) && (std::fabs(zi) <= tolODz))    // sd>0
                {
                  if (sd > dRmax) // Avoid rounding errors due to precision issue
                  {
                    // seen on 64 bits systems. Split and recompute
                    double fTerm = sd - std::fmod(sd, dRmax);
                    sd = fTerm + DistanceToIn(p + fTerm * v, v);
                  }
                  if (!fPhiFullCone)
                  {
                    xi     = p.x() + sd * v.x();
                    yi     = p.y() + sd * v.y();
                    ri     = rMinAv + zi * tanRMin;
                    cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

                    if (cosPsi >= cosHDPhiIT)
                    {
                      snxt = sd;
                    }
                  }
                  else
                  {
                    return sd;
                  }
                }
              }
              else
              {
                return UUtils::kInfinity;
              }
            }
          }
        }
        else   // 2nd root
        {
          b = nt2 / nt1;
          c = nt3 / nt1;
          d = b * b - c;

          if (d > 0)
          {
            if (b > 0)
            {
              sd = c / (-b - std::sqrt(d));
            }
            else
            {
              sd = -b + std::sqrt(d);
            }
            zi = p.z() + sd * v.z();

            if ((sd >= 0) && (std::fabs(zi) <= tolODz))    // sd>0
            {
              if (sd > dRmax) // Avoid rounding errors due to precision issues
              {
                // seen on 64 bits systems. Split and recompute
                double fTerm = sd - std::fmod(sd, dRmax);
                sd = fTerm + DistanceToIn(p + fTerm * v, v);
              }
              if (!fPhiFullCone)
              {
                xi     = p.x() + sd * v.x();
                yi     = p.y() + sd * v.y();
                ri     = rMinAv + zi * tanRMin;
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / ri;

                if (cosPsi >= cosHDPhiIT)
                {
                  snxt = sd;
                }
              }
              else
              {
                return sd;
              }
            }
          }
        }
      }
    }
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to VUSolid::Tolerance()*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> Should use some form of loop Construct

  if (!fPhiFullCone)
  {
    // First phi surface (starting phi)

    Comp    = v.x() * sinSPhi - v.y() * cosSPhi;

    if (Comp < 0)      // Component in outwards normal dirn
    {
      Dist = (p.y() * cosSPhi - p.x() * sinSPhi);

      if (Dist < halfCarTolerance)
      {
        sd = Dist / Comp;

        if (sd < snxt)
        {
          if (sd < 0)
          {
            sd = 0.0;
          }

          zi = p.z() + sd * v.z();

          if (std::fabs(zi) <= tolODz)
          {
            xi        = p.x() + sd * v.x();
            yi        = p.y() + sd * v.y();
            rhoi2    = xi * xi + yi * yi;
            tolORMin2 = (rMinOAv + zi * tanRMin) * (rMinOAv + zi * tanRMin);
            tolORMax2 = (rMaxOAv + zi * tanRMax) * (rMaxOAv + zi * tanRMax);

            if ((rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2))
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ((yi * cosCPhi - xi * sinCPhi) <= 0)
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }

    // Second phi surface (Ending phi)

    Comp    = -(v.x() * sinEPhi - v.y() * cosEPhi);

    if (Comp < 0)     // Component in outwards normal dirn
    {
      Dist = -(p.y() * cosEPhi - p.x() * sinEPhi);
      if (Dist < halfCarTolerance)
      {
        sd = Dist / Comp;

        if (sd < snxt)
        {
          if (sd < 0)
          {
            sd = 0.0;
          }

          zi = p.z() + sd * v.z();

          if (std::fabs(zi) <= tolODz)
          {
            xi        = p.x() + sd * v.x();
            yi        = p.y() + sd * v.y();
            rhoi2    = xi * xi + yi * yi;
            tolORMin2 = (rMinOAv + zi * tanRMin) * (rMinOAv + zi * tanRMin);
            tolORMax2 = (rMaxOAv + zi * tanRMax) * (rMaxOAv + zi * tanRMax);

            if ((rhoi2 >= tolORMin2) && (rhoi2 <= tolORMax2))
            {
              // z and r intersections good - check intersecting with
              // correct half-plane

              if ((yi * cosCPhi - xi * sinCPhi) >= 0.0)
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }
  }
  if (snxt < halfCarTolerance)
  {
    snxt = 0.;
  }

  return snxt;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

double UCons::SafetyFromOutside(const UVector3& p, bool aAccurate) const
{
  double safe = 0.0, rho, safeR1, safeR2, safeZ, safePhi;
  double pRMin, pRMax;
  bool outside;

  rho  = std::sqrt(p.x() * p.x() + p.y() * p.y());
  safeZ = std::fabs(p.z()) - fDz;
  safeR1 = 0; safeR2 = 0;

  if (fRmin1 || fRmin2)
  {
    pRMin  = tanRMin * p.z() + (fRmin1 + fRmin2) * 0.5;
    safeR1  = (pRMin - rho) / secRMin;

    pRMax  = tanRMax * p.z() + (fRmax1 + fRmax2) * 0.5;
    safeR2  = (rho - pRMax) / secRMax;

    if (safeR1 > safeR2)
    {
      safe = safeR1;
    }
    else
    {
      safe = safeR2;
    }
  }
  else
  {
    pRMax  = tanRMax * p.z() + (fRmax1 + fRmax2) * 0.5;
    safe    = (rho - pRMax) / secRMax;
  }
  if (safeZ > safe)
  {
    safe = safeZ;
  }

  if (!fPhiFullCone && rho)
  {
    // Psi=angle from central phi to point
    //
    safePhi = SafetyToPhi(p,rho,outside);
    if ((outside) && (safePhi > safe))
    {
      safe = safePhi;
    }
  }
  if (safe < 0.0)
  {
    safe = 0.0; return safe; //point is Inside
  }
  if (!aAccurate) return safe;

  double safsq = 0.0;
  int count = 0;
  if (safeR1 > 0)
  {
    safsq += safeR1 * safeR1;
    count++;
  }
  if (safeR2 > 0)
  {
    safsq += safeR2 * safeR2;
    count++;
  }
  if (safeZ > 0)
  {
    safsq += safeZ * safeZ;
    count++;
  }
  if (count == 1) return safe;
  return std::sqrt(safsq);
}

///////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection




// double DistanceToOut( const UVector3& p, const UVector3& v, bool calcNorm, bool aConvex, UVector3 *n ) const;

double UCons::DistanceToOut(const UVector3& p,
                            const UVector3&  v,
                            UVector3&       aNormalVector,
                            bool&           aConvex,
                            double /* aPstep*/) const
{
  ESide side = kNull, sider = kNull, sidephi = kNull;

  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfRadTolerance = kRadTolerance * 0.5;
  static const double halfAngTolerance = kAngTolerance * 0.5;

  double snxt, srd, sphi, pdist;

  double rMaxAv;  // Data for outer cone
  double rMinAv;  // Data for inner cone

  double t1, t2, t3, rout, rin, nt1, nt2, nt3;
  double b, c, d, sr2, sr3;

  // Vars for intersection within tolerance

  ESide   sidetol = kNull;
  double slentol = UUtils::kInfinity;

  // Vars for phi intersection:

  double pDistS, compS, pDistE, compE, sphi2, xi, yi, risec, vphi;
  double zi, ri, deltaRoi2;

  // Z plane intersection

  if (v.z() > 0.0)
  {
    pdist = fDz - p.z();

    if (pdist > halfCarTolerance)
    {
      snxt = pdist / v.z();
      side = kPZ;
    }
    else
    {
      aNormalVector        = UVector3(0, 0, 1);
      aConvex = true;
      return  snxt = 0.0;
    }
  }
  else if (v.z() < 0.0)
  {
    pdist = fDz + p.z();

    if (pdist > halfCarTolerance)
    {
      snxt = -pdist / v.z();
      side = kMZ;
    }
    else
    {
      aNormalVector        = UVector3(0, 0, -1);
      aConvex = true;
      return snxt = 0.0;
    }
  }
  else     // Travel perpendicular to z axis
  {
    snxt = UUtils::kInfinity;
    side = kNull;
  }

  // Radial Intersections
  //
  // Intersection with outer cone (possible return) and
  //                   inner cone (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x()+t*v.x() etc.
  //
  // Intersects with x^2+y^2=(a*z+b)^2
  //
  // where a=tanRMax or tanRMin
  //       b=rMaxAv or rMinAv
  //
  // (vx^2+vy^2-(a*vz)^2)t^2+2t(pxvx+pyvy-a*vz(a*pz+b))+px^2+py^2-(a*pz+b)^2=0;
  //     t1                       t2                      t3
  //
  //  \--------u-------/       \-----------v----------/ \---------w--------/

  rMaxAv  = (fRmax1 + fRmax2) * 0.5;

  t1   = 1.0 - v.z() * v.z();   // since v normalised
  t2   = p.x() * v.x() + p.y() * v.y();
  t3   = p.x() * p.x() + p.y() * p.y();
  rout = tanRMax * p.z() + rMaxAv;

  nt1 = t1 - (tanRMax * v.z()) * (tanRMax * v.z());
  nt2 = t2 - tanRMax * v.z() * rout;
  nt3 = t3 - rout * rout;

  if (v.z() > 0.0)
  {
    deltaRoi2 = snxt * snxt * t1 + 2 * snxt * t2 + t3
                - fRmax2 * (fRmax2 + kRadTolerance * secRMax);
  }
  else if (v.z() < 0.0)
  {
    deltaRoi2 = snxt * snxt * t1 + 2 * snxt * t2 + t3
                - fRmax1 * (fRmax1 + kRadTolerance * secRMax);
  }
  else
  {
    deltaRoi2 = 1.0;
  }

  if (nt1 && (deltaRoi2 > 0.0))
  {
    // Equation quadratic => 2 roots : second root must be leaving

    b = nt2 / nt1;
    c = nt3 / nt1;
    d = b * b - c;

    if (d >= 0)
    {
      // Check if on outer cone & heading outwards
      // NOTE: Should use rho-rout>-kRadTolerance*0.5

      if (nt3 > -halfRadTolerance && nt2 >= 0)
      {
        risec     = std::sqrt(t3) * secRMax;
        aConvex = true;
        aNormalVector        = UVector3(p.x() / risec, p.y() / risec, -tanRMax / secRMax);
        return snxt = 0;
      }
      else
      {
        sider = kRMax ;
        if (b > 0)
        {
          srd = -b - std::sqrt(d);
        }
        else
        {
          srd = c / (-b + std::sqrt(d));
        }

        zi    = p.z() + srd * v.z();
        ri    = tanRMax * zi + rMaxAv;

        if ((ri >= 0) && (-halfRadTolerance <= srd) && (srd <= halfRadTolerance))
        {
          // An intersection within the tolerance
          //   we will Store it in case it is good -
          //
          slentol = srd;
          sidetol = kRMax;
        }
        if ((ri < 0) || (srd < halfRadTolerance))
        {
          // Safety: if both roots -ve ensure that srd cannot `win'
          //         distance to out

          if (b > 0)
          {
            sr2 = c / (-b - std::sqrt(d));
          }
          else
          {
            sr2 = -b + std::sqrt(d);
          }
          zi  = p.z() + sr2 * v.z();
          ri  = tanRMax * zi + rMaxAv;

          if ((ri >= 0) && (sr2 > halfRadTolerance))
          {
            srd = sr2;
          }
          else
          {
            srd = UUtils::kInfinity;

            if ((-halfRadTolerance <= sr2) && (sr2 <= halfRadTolerance))
            {
              // An intersection within the tolerance.
              // Storing it in case it is good.

              slentol = sr2;
              sidetol = kRMax;
            }
          }
        }
      }
    }
    else
    {
      // No intersection with outer cone & not parallel
      // -> already outside, no intersection

      risec     = std::sqrt(t3) * secRMax;
      aConvex = true;
      aNormalVector        = UVector3(p.x() / risec, p.y() / risec, -tanRMax / secRMax);
      return snxt = 0.0;
    }
  }
  else if (nt2 && (deltaRoi2 > 0.0))
  {
    // Linear case (only one intersection) => point outside outer cone

    risec     = std::sqrt(t3) * secRMax;
    aConvex = true;
    aNormalVector        = UVector3(p.x() / risec, p.y() / risec, -tanRMax / secRMax);
    return snxt = 0.0;
  }
  else
  {
    // No intersection -> parallel to outer cone
    // => Z or inner cone intersection

    srd = UUtils::kInfinity;
  }

  // Check possible intersection within tolerance

  if (slentol <= halfCarTolerance)
  {
    // An intersection within the tolerance was found.
    // We must accept it only if the momentum points outwards.
    //
    // UVector3 ptTol;  // The point of the intersection
    // ptTol= p + slentol*v;
    // ri=tanRMax*zi+rMaxAv;
    //
    // Calculate a normal vector, as below

    xi    = p.x() + slentol * v.x();
    yi    = p.y() + slentol * v.y();
    risec = std::sqrt(xi * xi + yi * yi) * secRMax;
    UVector3 norm = UVector3(xi / risec, yi / risec, -tanRMax / secRMax);

    if (norm.Dot(v) > 0)     // We will leave the Cone immediatelly
    {
      aNormalVector        = norm.Unit();
      aConvex = true;
      return snxt = 0.0;
    }
    else // On the surface, but not heading out so we ignore this intersection
    {
      //                                        (as it is within tolerance).
      slentol = UUtils::kInfinity;
    }
  }

  // Inner Cone intersection

  if (fRmin1 || fRmin2)
  {
    nt1    = t1 - (tanRMin * v.z()) * (tanRMin * v.z());

    if (nt1)
    {
      rMinAv  = (fRmin1 + fRmin2) * 0.5;
      rin    = tanRMin * p.z() + rMinAv;
      nt2    = t2 - tanRMin * v.z() * rin;
      nt3    = t3 - rin * rin;

      // Equation quadratic => 2 roots : first root must be leaving

      b = nt2 / nt1;
      c = nt3 / nt1;
      d = b * b - c;

      if (d >= 0.0)
      {
        // NOTE: should be rho-rin<kRadTolerance*0.5,
        //       but using squared versions for efficiency

        if (nt3 < kRadTolerance * (rin + kRadTolerance * 0.25))
        {
          if (nt2 < 0.0)
          {
            aConvex = false;
            risec = std::sqrt(p.x() * p.x() + p.y() * p.y()) * secRMin;
            aNormalVector = UVector3(-p.x() / risec, -p.y() / risec, tanRMin / secRMin);
            return          snxt      = 0.0;
          }
        }
        else
        {
          if (b > 0)
          {
            sr2 = -b - std::sqrt(d);
          }
          else
          {
            sr2 = c / (-b + std::sqrt(d));
          }
          zi  = p.z() + sr2 * v.z();
          ri  = tanRMin * zi + rMinAv;

          if ((ri >= 0.0) && (-halfRadTolerance <= sr2) && (sr2 <= halfRadTolerance))
          {
            // An intersection within the tolerance
            // storing it in case it is good.

            slentol = sr2;
            sidetol = kRMax;
          }
          if ((ri < 0) || (sr2 < halfRadTolerance))
          {
            if (b > 0)
            {
              sr3 = c / (-b - std::sqrt(d));
            }
            else
            {
              sr3 = -b + std::sqrt(d);
            }

            // Safety: if both roots -ve ensure that srd cannot `win'
            //         distancetoout

            if (sr3 > halfRadTolerance)
            {
              if (sr3 < srd)
              {
                zi = p.z() + sr3 * v.z();
                ri = tanRMin * zi + rMinAv;

                if (ri >= 0.0)
                {
                  srd = sr3;
                  sider = kRMin;
                }
              }
            }
            else if (sr3 > -halfRadTolerance)
            {
              // Intersection in tolerance. Store to check if it's good

              slentol = sr3;
              sidetol = kRMin;
            }
          }
          else if ((sr2 < srd) && (sr2 > halfCarTolerance))
          {
            srd  = sr2;
            sider = kRMin;
          }
          else if (sr2 > -halfCarTolerance)
          {
            // Intersection in tolerance. Store to check if it's good

            slentol = sr2;
            sidetol = kRMin;
          }
          if (slentol <= halfCarTolerance)
          {
            // An intersection within the tolerance was found.
            // We must accept it only if  the momentum points outwards.

            UVector3 norm;

            // Calculate a normal vector, as below

            xi     = p.x() + slentol * v.x();
            yi     = p.y() + slentol * v.y();
            if (sidetol == kRMax)
            {
              risec = std::sqrt(xi * xi + yi * yi) * secRMax;
              norm = UVector3(xi / risec, yi / risec, -tanRMax / secRMax);
            }
            else
            {
              risec = std::sqrt(xi * xi + yi * yi) * secRMin;
              norm = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
            }
            if (norm.Dot(v) > 0)
            {
              // We will leave the cone immediately

              aNormalVector        = norm.Unit();
              aConvex = true;
              return snxt = 0.0;
            }
            else
            {
              // On the surface, but not heading out so we ignore this
              // intersection (as it is within tolerance).

              slentol = UUtils::kInfinity;
            }
          }
        }
      }
    }
  }

  // Linear case => point outside inner cone ---> outer cone intersect
  //
  // Phi Intersection

  if (!fPhiFullCone)
  {
    // add angle calculation with correction
    // of the difference in domain of atan2 and Sphi

    vphi = std::atan2(v.y(), v.x());

    if (vphi < fSPhi - halfAngTolerance)
    {
      vphi += 2 * UUtils::kPi;
    }
    else if (vphi > fSPhi + fDPhi + halfAngTolerance)
    {
      vphi -= 2 * UUtils::kPi;
    }

    if (p.x() || p.y())     // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS = p.x() * sinSPhi - p.y() * cosSPhi;
      pDistE = -p.x() * sinEPhi + p.y() * cosEPhi;

      // Comp -ve when in direction of outwards normal

      compS = -sinSPhi * v.x() + cosSPhi * v.y();
      compE = sinEPhi * v.x() - cosEPhi * v.y();

      sidephi = kNull;

      if (((fDPhi <= UUtils::kPi) && ((pDistS <= halfCarTolerance)
                                      && (pDistE <= halfCarTolerance)))
          || ((fDPhi > UUtils::kPi) && !((pDistS > halfCarTolerance)
                                         && (pDistE >  halfCarTolerance))))
      {
        // Inside both phi *full* planes
        if (compS < 0)
        {
          sphi = pDistS / compS;
          if (sphi >= -halfCarTolerance)
          {
            xi = p.x() + sphi * v.x();
            yi = p.y() + sphi * v.y();

            // Check intersecting with correct half-plane
            // (if not -> no intersect)
            //
            if ((std::fabs(xi) <= VUSolid::Tolerance())
                && (std::fabs(yi) <= VUSolid::Tolerance()))
            {
              sidephi = kSPhi;
              if ((fSPhi - halfAngTolerance <= vphi)
                  && (fSPhi + fDPhi + halfAngTolerance >= vphi))
              {
                sphi = UUtils::kInfinity;
              }
            }
            else if ((yi * cosCPhi - xi * sinCPhi) >= 0)
            {
              sphi = UUtils::kInfinity;
            }
            else
            {
              sidephi = kSPhi;
              if (pDistS > -halfCarTolerance)
              {
                sphi = 0.0; // Leave by sphi immediately
              }
            }
          }
          else
          {
            sphi = UUtils::kInfinity;
          }
        }
        else
        {
          sphi = UUtils::kInfinity;
        }

        if (compE < 0)
        {
          sphi2 = pDistE / compE;

          // Only check further if < starting phi intersection
          //
          if ((sphi2 > -halfCarTolerance) && (sphi2 < sphi))
          {
            xi = p.x() + sphi2 * v.x();
            yi = p.y() + sphi2 * v.y();

            // Check intersecting with correct half-plane

            if ((std::fabs(xi) <= VUSolid::Tolerance())
                && (std::fabs(yi) <= VUSolid::Tolerance()))
            {
              // Leaving via ending phi

              if (!((fSPhi - halfAngTolerance <= vphi)
                    && (fSPhi + fDPhi + halfAngTolerance >= vphi)))
              {
                sidephi = kEPhi;
                if (pDistE <= -halfCarTolerance)
                {
                  sphi = sphi2;
                }
                else
                {
                  sphi = 0.0;
                }
              }
            }
            else // Check intersecting with correct half-plane
              if (yi * cosCPhi - xi * sinCPhi >= 0)
              {
                // Leaving via ending phi

                sidephi = kEPhi;
                if (pDistE <= -halfCarTolerance)
                {
                  sphi = sphi2;
                }
                else
                {
                  sphi = 0.0;
                }
              }
          }
        }
      }
      else
      {
        sphi = UUtils::kInfinity;
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if ((fSPhi - halfAngTolerance <= vphi)
          && (vphi <= fSPhi + fDPhi + halfAngTolerance))
      {
        sphi = UUtils::kInfinity;
      }
      else
      {
        sidephi = kSPhi ;  // arbitrary
        sphi    = 0.0;
      }
    }
    if (sphi < snxt)   // Order intersecttions
    {
      snxt = sphi;
      side = sidephi;
    }
  }
  if (srd < snxt)    // Order intersections
  {
    snxt = srd  ;
    side = sider;
  }

  switch (side)
  {
      // Note: returned vector not normalised
    case kRMax:        // (divide by frmax for Unit vector)
      xi         = p.x() + snxt * v.x();
      yi         = p.y() + snxt * v.y();
      risec     = std::sqrt(xi * xi + yi * yi) * secRMax;
      aNormalVector        = UVector3(xi / risec, yi / risec, -tanRMax / secRMax);
      aConvex = true;
      break;
    case kRMin:
      xi         = p.x() + snxt * v.x();
      yi         = p.y() + snxt * v.y();
      risec = std::sqrt(xi * xi + yi * yi) * secRMin;
      aNormalVector = UVector3(-xi / risec, -yi / risec, tanRMin / secRMin);
      aConvex = false;  // Rmin is inconvex
      break;
    case kSPhi:
      if (fDPhi <= UUtils::kPi)
      {
        aNormalVector        = UVector3(sinSPhi, -cosSPhi, 0);
        aConvex = true;
      }
      else
      {
        aNormalVector = UVector3(sinSPhi, -cosSPhi, 0);
        aConvex = false;
      }
      break;
    case kEPhi:
      if (fDPhi <= UUtils::kPi)
      {
        aNormalVector = UVector3(-sinEPhi, cosEPhi, 0);
        aConvex = true;
      }
      else
      {
        aNormalVector = UVector3(-sinEPhi, cosEPhi, 0);
        aConvex = false;
      }
      break;
    case kPZ:
      aNormalVector        = UVector3(0, 0, 1);
      aConvex = true;
      break;
    case kMZ:
      aNormalVector        = UVector3(0, 0, -1);
      aConvex = true;
      break;
    default:
      cout << std::endl;
      // DumpInfo();
      std::ostringstream message;
      int oldprc = message.precision(16);
      message << "Undefined side for valid surface normal to solid."
              << std::endl
              << "Position:"  << std::endl << std::endl
              << "p.x = "  << p.x() << " mm" << std::endl
              << "p.y = "  << p.y() << " mm" << std::endl
              << "p.z = "  << p.z() << " mm" << std::endl << std::endl
              << "pho at z = "   << std::sqrt(p.x() * p.x() + p.y() * p.y())
              << " mm" << std::endl << std::endl;
      if (p.x() != 0. || p.y() != 0.)
      {
        message << "point phi = "   << std::atan2(p.y(), p.x()) / (UUtils::kPi / 180.0)
                << " degree" << std::endl << std::endl;
      }
      message << "Direction:" << std::endl << std::endl
              << "v.x = "  << v.x() << std::endl
              << "v.y = "  << v.y() << std::endl
              << "v.z = "  << v.z() << std::endl << std::endl
              << "Proposed distance :" << std::endl << std::endl
              << "snxt = "    << snxt << " mm" << std::endl;
      message.precision(oldprc);
      UUtils::Exception("UCons::DistanceToOut()", "GeomSolids0002", UWarning, 1, message.str().c_str());
      break;
  }
  if (snxt < halfCarTolerance)
  {
    snxt = 0.;
  }

  return snxt;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

double UCons::SafetyFromInside(const UVector3& p, bool) const
{
  double safe = 0.0, rho, safeZ;
 
#ifdef UCSGDEBUG
  if (Inside(p) == eOutside)
  {
    int oldprc = cout.precision(16);
    cout << std::endl;
    DumpInfo();
    cout << "Position:" << std::endl << std::endl;
    cout << "p.x = "   << p.x() << " mm" << std::endl;
    cout << "p.y = "   << p.y() << " mm" << std::endl;
    cout << "p.z = "   << p.z() << " mm" << std::endl << std::endl;
    cout << "pho at z = "  << std::sqrt(p.x() * p.x() + p.y() * p.y())
         << " mm" << std::endl << std::endl;
    if ((p.x() != 0.) || (p.x() != 0.))
    {
      cout << "point phi = "   << std::atan2(p.y(), p.x()) / degree
           << " degree" << std::endl << std::endl;
    }
    cout.precision(oldprc);
    UUtils::Exception("UCons::UCons()", "GeomSolids0002", UWarning, 1, message.str().c_str());

  }
#endif
  
  rho = std::sqrt(p.x() * p.x() + p.y() * p.y());

  safe = SafetyFromInsideR(p, rho, false);
  safeZ = fDz - std::fabs(p.z());

  if (safeZ < safe)
  {
    safe = safeZ;
  }
  if (safe < 0)
  {
    safe = 0;
  }

  return safe;
}


//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

UGeometryType UCons::GetEntityType() const
{
  return std::string("Cons");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
VUSolid* UCons::Clone() const
{

  return new UCons(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& UCons::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: UCons\n"
     << " Parameters: \n"
     << "         inside -fDz radius : " << fRmin1 << " mm \n"
     << "         outside -fDz radius: " << fRmax1 << " mm \n"
     << "         inside +fDz radius : " << fRmin2 << " mm \n"
     << "         outside +fDz radius: " << fRmax2 << " mm \n"
     << "         half length in Z   : " << fDz << " mm \n"
     << "         starting angle of segment: " << fSPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "         delta angle of segment   : " << fDPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}



/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

UVector3 UCons::GetPointOnSurface() const
{
  // declare working variables
  //
  double Aone, Atwo, Athree, Afour, Afive, slin, slout, phi;
  double zRand, cosu, sinu, rRand1, rRand2, chose, rone, rtwo, qone, qtwo;
  rone = (fRmax1 - fRmax2) / (2.*fDz);
  rtwo = (fRmin1 - fRmin2) / (2.*fDz);
  qone = 0.;
  qtwo = 0.;
  if (fRmax1 != fRmax2)
  {
    qone = fDz * (fRmax1 + fRmax2) / (fRmax1 - fRmax2);
  }
  if (fRmin1 != fRmin2)
  {
    qtwo = fDz * (fRmin1 + fRmin2) / (fRmin1 - fRmin2);
  }
  slin   = std::sqrt(UUtils::sqr(fRmin1 - fRmin2) + UUtils::sqr(2.*fDz));
  slout = std::sqrt(UUtils::sqr(fRmax1 - fRmax2) + UUtils::sqr(2.*fDz));
  Aone   = 0.5 * fDPhi * (fRmax2 + fRmax1) * slout;
  Atwo   = 0.5 * fDPhi * (fRmin2 + fRmin1) * slin;
  Athree = 0.5 * fDPhi * (fRmax1 * fRmax1 - fRmin1 * fRmin1);
  Afour = 0.5 * fDPhi * (fRmax2 * fRmax2 - fRmin2 * fRmin2);
  Afive = fDz * (fRmax1 - fRmin1 + fRmax2 - fRmin2);

  phi   = UUtils::Random(fSPhi, fSPhi + fDPhi);
  cosu   = std::cos(phi);
  sinu = std::sin(phi);
  rRand1 = UUtils::GetRadiusInRing(fRmin1, fRmax1);
  rRand2 = UUtils::GetRadiusInRing(fRmin2, fRmax2);
 
  if ((fSPhi == 0.) && fPhiFullCone)
  {
    Afive = 0.;
  }
  chose = UUtils::Random(0., Aone + Atwo + Athree + Afour + 2.*Afive);

  if ((chose >= 0.) && (chose < Aone))
  {
    if (fRmin1 != fRmin2)
    {
      zRand = UUtils::Random(-1.*fDz, fDz);
      return UVector3(rtwo * cosu * (qtwo - zRand),
                      rtwo * sinu * (qtwo - zRand), zRand);
    }
    else
    {
      return UVector3(fRmin1 * cosu, fRmin2 * sinu,
                      UUtils::Random(-1.*fDz, fDz));
    }
  }
  else if ((chose >= Aone) && (chose <= Aone + Atwo))
  {
    if (fRmax1 != fRmax2)
    {
      zRand = UUtils::Random(-1.*fDz, fDz);
      return UVector3(rone * cosu * (qone - zRand),
                      rone * sinu * (qone - zRand), zRand);
    }
    else
    {
      return UVector3(fRmax1 * cosu, fRmax2 * sinu,
                      UUtils::Random(-1.*fDz, fDz));
    }
  }
  else if ((chose >= Aone + Atwo) && (chose < Aone + Atwo + Athree))
  {
    return UVector3(rRand1 * cosu, rRand1 * sinu, -1 * fDz);
  }
  else if ((chose >= Aone + Atwo + Athree)
           && (chose < Aone + Atwo + Athree + Afour))
  {
    return UVector3(rRand2 * cosu, rRand2 * sinu, fDz);
  }
  else if ((chose >= Aone + Atwo + Athree + Afour)
           && (chose < Aone + Atwo + Athree + Afour + Afive))
  {
    zRand = UUtils::Random(-1.*fDz, fDz);
    rRand1 = UUtils::Random(fRmin2 - ((zRand - fDz) / (2.*fDz)) * (fRmin1 - fRmin2),
                            fRmax2 - ((zRand - fDz) / (2.*fDz)) * (fRmax1 - fRmax2));
    return UVector3(rRand1 * std::cos(fSPhi),
                    rRand1 * std::sin(fSPhi), zRand);
  }
  else
  {
    zRand = UUtils::Random(-1.*fDz, fDz);
    rRand1 = UUtils::Random(fRmin2 - ((zRand - fDz) / (2.*fDz)) * (fRmin1 - fRmin2),
                            fRmax2 - ((zRand - fDz) / (2.*fDz)) * (fRmax1 - fRmax2));
    return UVector3(rRand1 * std::cos(fSPhi + fDPhi),
                    rRand1 * std::sin(fSPhi + fDPhi), zRand);
  }
}


void UCons::Extent(UVector3& aMin, UVector3& aMax) const
{
  double max = fRmax1 > fRmax2 ? fRmax1 : fRmax2;
  aMin = UVector3(-max, -max, -fDz);
  aMax = UVector3(max, max, fDz);
}

void UCons::CheckDPhiAngle(double dPhi)
{
  fPhiFullCone = true;
  if (dPhi >= 2 * UUtils::kPi - kAngTolerance * 0.5)
  {
    fDPhi = 2 * UUtils::kPi;
    fSPhi = 0;
  }
  else
  {
    fPhiFullCone = false;
    if (dPhi > 0)
    {
      fDPhi = dPhi;
    }
    else
    {
      std::ostringstream message;
      message << "Invalid dphi." << std::endl
              << "Negative or zero delta-Phi (" << dPhi << ") in solid: "
              << GetName();
      UUtils::Exception("UCons::CheckDPhiAngle()", "GeomSolids0002",
                        UFatalErrorInArguments, 1, message.str().c_str());
    }
  }
}

void UCons::GetParametersList(int, double* aArray) const
{
  aArray[0] = GetInnerRadiusMinusZ();
  aArray[1] = GetOuterRadiusMinusZ();
  aArray[2] = GetInnerRadiusPlusZ();
  aArray[3] = GetOuterRadiusPlusZ();
  aArray[4] = GetZHalfLength();
  aArray[5] = GetStartPhiAngle();
  aArray[6] = GetDeltaPhiAngle();
}
