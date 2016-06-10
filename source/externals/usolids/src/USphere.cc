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
// USphere
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "USphere.hh"

using namespace std;

// Private enum: Not for external use - used by distanceToOut

enum ESide {kNull, kRMin, kRMax, kSPhi, kEPhi, kSTheta, kETheta};

// used by normal

enum ENorm {kNRMin, kNRMax, kNSPhi, kNEPhi, kNSTheta, kNETheta};

////////////////////////////////////////////////////////////////////////
//
// constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pDPhi>2PI then reset to 2PI

USphere::USphere(const std::string& pName,
                 double pRmin, double pRmax,
                 double pSPhi, double pDPhi,
                 double pSTheta, double pDTheta)
  : VUSolid(pName), fCubicVolume(0.),
    fSurfaceArea(0.),fEpsilon(2.e-11),
    fFullPhiSphere(true), fFullThetaSphere(true)
{
  kAngTolerance = faTolerance;

  // Check radii and Set radial tolerances

  kRadTolerance = frTolerance;
  if ((pRmin >= pRmax) || (pRmax < 1.1 * kRadTolerance) || (pRmin < 0))
  {
    std::ostringstream message;
    message << "Invalid radii for Solid: " << GetName() << std::endl
            << "pRmin = " << pRmin << ", pRmax = " << pRmax;
    UUtils::Exception("USphere::USphere()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fRmin = pRmin;
  fRmax = pRmax;
  fRminTolerance = (fRmin) ? std::max(kRadTolerance, fEpsilon * fRmin) : 0;
  kTolerance = std::max(kRadTolerance, fEpsilon * fRmax);

  // Check angles

  CheckPhiAngles(pSPhi, pDPhi);
  CheckThetaAngles(pSTheta, pDTheta);
}

/////////////////////////////////////////////////////////////////////
//
// Destructor

USphere::~USphere()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

USphere::USphere(const USphere& rhs)
  : VUSolid(rhs), fCubicVolume(0.),fSurfaceArea(0.),
    fRminTolerance(rhs.fRminTolerance),
    kTolerance(rhs.kTolerance), kAngTolerance(rhs.kAngTolerance),
    kRadTolerance(rhs.kRadTolerance), fEpsilon(rhs.fEpsilon),
    fRmin(rhs.fRmin), fRmax(rhs.fRmax), fSPhi(rhs.fSPhi), fDPhi(rhs.fDPhi),
    fSTheta(rhs.fSTheta), fDTheta(rhs.fDTheta),
    sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi),
    sinEPhi(rhs.sinEPhi), cosEPhi(rhs.cosEPhi),
    hDPhi(rhs.hDPhi), cPhi(rhs.cPhi), ePhi(rhs.ePhi),
    sinSTheta(rhs.sinSTheta), cosSTheta(rhs.cosSTheta),
    sinETheta(rhs.sinETheta), cosETheta(rhs.cosETheta),
    tanSTheta(rhs.tanSTheta), tanSTheta2(rhs.tanSTheta2),
    tanETheta(rhs.tanETheta), tanETheta2(rhs.tanETheta2), eTheta(rhs.eTheta),
    fFullPhiSphere(rhs.fFullPhiSphere), fFullThetaSphere(rhs.fFullThetaSphere),
    fFullSphere(rhs.fFullSphere)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

USphere& USphere::operator = (const USphere& rhs)
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
  fCubicVolume = rhs.fCubicVolume;
  fSurfaceArea = rhs.fSurfaceArea;
  fRminTolerance = rhs.fRminTolerance;
  kTolerance = rhs.kTolerance;
  kAngTolerance = rhs.kAngTolerance;
  kRadTolerance = rhs.kRadTolerance;
  fEpsilon = rhs.fEpsilon;
  fRmin = rhs.fRmin;
  fRmax = rhs.fRmax;
  fSPhi = rhs.fSPhi;
  fDPhi = rhs.fDPhi;
  fSTheta = rhs.fSTheta;
  fDTheta = rhs.fDTheta;
  sinCPhi = rhs.sinCPhi;
  cosCPhi = rhs.cosCPhi;
  cosHDPhiOT = rhs.cosHDPhiOT;
  cosHDPhiIT = rhs.cosHDPhiIT;
  sinSPhi = rhs.sinSPhi;
  cosSPhi = rhs.cosSPhi;
  sinEPhi = rhs.sinEPhi;
  cosEPhi = rhs.cosEPhi;
  hDPhi = rhs.hDPhi;
  cPhi = rhs.cPhi;
  ePhi = rhs.ePhi;
  sinSTheta = rhs.sinSTheta;
  cosSTheta = rhs.cosSTheta;
  sinETheta = rhs.sinETheta;
  cosETheta = rhs.cosETheta;
  tanSTheta = rhs.tanSTheta;
  tanSTheta2 = rhs.tanSTheta2;
  tanETheta = rhs.tanETheta;
  tanETheta2 = rhs.tanETheta2;
  eTheta = rhs.eTheta;
  fFullPhiSphere = rhs.fFullPhiSphere;
  fFullThetaSphere = rhs.fFullThetaSphere;
  fFullSphere = rhs.fFullSphere;

  return *this;
}

///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
// Split into radius, phi, theta checks
// Each check modifies 'in', or returns as approprate

VUSolid::EnumInside USphere::Inside(const UVector3& p) const
{
  double rho, rho2, rad2, tolRMin, tolRMax;
  double pPhi, pTheta;
  VUSolid::EnumInside in = eOutside;
  static const double halfAngTolerance = kAngTolerance * 0.5;
  const double halfTolerance = kTolerance * 0.5;
  const double halfRminTolerance = fRminTolerance * 0.5;
  const double rMaxMinus = fRmax - halfTolerance;
  const double rMinPlus = (fRmin > 0) ? fRmin + halfRminTolerance : 0;

  rho2 = p.x() * p.x() + p.y() * p.y();
  rad2 = rho2 + p.z() * p.z();

  // Check radial surfaces. Sets 'in'

  tolRMin = rMinPlus;
  tolRMax = rMaxMinus;

  if(rad2 == 0.0)
  { 
    if (fRmin > 0.0)
    {
      return in = eOutside;
    }
    if ( (!fFullPhiSphere) || (!fFullThetaSphere) )
    {
      return in = eSurface;
    }
    else
    {
      return in = eInside; 
    }
  }

  if ((rad2 <= rMaxMinus * rMaxMinus) && (rad2 >= rMinPlus * rMinPlus))
  {
    in = eInside;
  }
  else
  {
    tolRMax = fRmax + halfTolerance;                  // outside case
    tolRMin = std::max(fRmin - halfRminTolerance, 0.);    // outside case
    if ((rad2 <= tolRMax * tolRMax) && (rad2 >= tolRMin * tolRMin))
    {
      in = eSurface;
    }
    else
    {
      return in = eOutside;
    }
  }

  // Phi boundaries  : Do not check if it has no phi boundary!

  if (!fFullPhiSphere && rho2)  // [fDPhi < 2*UUtils::kPi] and [p.x() or p.y()]
  {
    pPhi = std::atan2(p.y(), p.x());

    if (pPhi < fSPhi - halfAngTolerance)
    {
      pPhi += 2 * UUtils::kPi;
    }
    else if (pPhi > ePhi + halfAngTolerance)
    {
      pPhi -= 2 * UUtils::kPi;
    }

    if ((pPhi < fSPhi - halfAngTolerance)
        || (pPhi > ePhi + halfAngTolerance))
    {
      return in = eOutside;
    }

    else if (in == eInside) // else it's eSurface anyway already
    {
      if ((pPhi < fSPhi + halfAngTolerance)
          || (pPhi > ePhi - halfAngTolerance))
      {
        in = eSurface;
      }
    }
  }

  // Theta bondaries

  if ((rho2 || p.z()) && (!fFullThetaSphere))
  {
    rho   = std::sqrt(rho2);
    pTheta = std::atan2(rho, p.z());

    if (in == eInside)
    {
      if ( ((fSTheta > 0.0) && (pTheta < fSTheta + halfAngTolerance))
        || ((eTheta < UUtils::kPi) && (pTheta > eTheta - halfAngTolerance)) ) 
      {

       if ( (( (fSTheta>0.0)&&(pTheta>=fSTheta-halfAngTolerance) )
            || (fSTheta == 0.0) )
          && ((eTheta==UUtils::kPi)||(pTheta <= eTheta + halfAngTolerance) ) )
        {
          in = eSurface;
        }
        else
        {
          in = eOutside;
        }
      }
    }
    else
    {
      if ( ((fSTheta > 0.0)&&(pTheta < fSTheta - halfAngTolerance))
	   ||((eTheta < UUtils::kPi  )&&(pTheta > eTheta + halfAngTolerance)) )
      {
        in = eOutside;
      }
    }
  }
  return in;
}

/////////////////////////////////////////////////////////////////////
//
// Return Unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

bool USphere::Normal(const UVector3& p, UVector3& n) const
{
  int noSurfaces = 0;
  double rho, rho2, radius, pTheta, pPhi = 0.;
  double distRMin = UUtils::Infinity();
  double distSPhi = UUtils::Infinity(), distEPhi = UUtils::Infinity();
  double distSTheta = UUtils::Infinity(), distETheta = UUtils::Infinity();
  UVector3 nR, nPs, nPe, nTs, nTe, nZ(0., 0., 1.);
  UVector3 norm, sumnorm(0., 0., 0.);

  static const double halfCarTolerance = 0.5 * VUSolid::Tolerance();
  static const double halfAngTolerance = 0.5 * kAngTolerance;

  rho2 = p.x() * p.x() + p.y() * p.y();
  radius = std::sqrt(rho2 + p.z() * p.z());
  rho = std::sqrt(rho2);

  double    distRMax = std::fabs(radius - fRmax);
  if (fRmin)  distRMin = std::fabs(radius - fRmin);

  if (rho && !fFullSphere)
  {
    pPhi = std::atan2(p.y(), p.x());

    if (pPhi < fSPhi - halfAngTolerance)
    {
      pPhi += 2 * UUtils::kPi;
    }
    else if (pPhi > ePhi + halfAngTolerance)
    {
      pPhi -= 2 * UUtils::kPi;
    }
  }
  if (!fFullPhiSphere)
  {
    if (rho)
    {
      distSPhi = std::fabs(pPhi - fSPhi);
      distEPhi = std::fabs(pPhi - ePhi);
    }
    else if (!fRmin)
    {
      distSPhi = 0.;
      distEPhi = 0.;
    }
    nPs = UVector3(sinSPhi, -cosSPhi, 0);
    nPe = UVector3(-sinEPhi, cosEPhi, 0);
  }
  if (!fFullThetaSphere)
  {
    if (rho)
    {
      pTheta     = std::atan2(rho, p.z());
      distSTheta = std::fabs(pTheta - fSTheta);
      distETheta = std::fabs(pTheta - eTheta);

      nTs = UVector3(-cosSTheta * p.x() / rho,
                     -cosSTheta * p.y() / rho,
                     sinSTheta);

      nTe = UVector3(cosETheta * p.x() / rho,
                     cosETheta * p.y() / rho,
                     -sinETheta);
    }
    else if (!fRmin)
    {
      if (fSTheta)
      {
        distSTheta = 0.;
        nTs = UVector3(0., 0., -1.);
      }
      if (eTheta < UUtils::kPi)
      {
        distETheta = 0.;
        nTe = UVector3(0., 0., 1.);
      }
    }
  }
  if (radius)
  {
    nR = UVector3(p.x() / radius, p.y() / radius, p.z() / radius);
  }

  if (distRMax <= halfCarTolerance)
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if (fRmin && (distRMin <= halfCarTolerance))
  {
    noSurfaces ++;
    sumnorm -= nR;
  }
  if (!fFullPhiSphere)
  {
    if (distSPhi <= halfAngTolerance)
    {
      noSurfaces ++;
      sumnorm += nPs;
    }
    if (distEPhi <= halfAngTolerance)
    {
      noSurfaces ++;
      sumnorm += nPe;
    }
  }
  if (!fFullThetaSphere)
  {
    if ((distSTheta <= halfAngTolerance) && (fSTheta > 0.))
    {
      noSurfaces ++;
      if ((radius <= halfCarTolerance) && fFullPhiSphere)
      {
        sumnorm += nZ;
      }
      else
      {
        sumnorm += nTs;
      }
    }
    if ((distETheta <= halfAngTolerance) && (eTheta < UUtils::kPi))
    {
      noSurfaces ++;
      if ((radius <= halfCarTolerance) && fFullPhiSphere)
      {
        sumnorm -= nZ;
      }
      else
      {
        sumnorm += nTe;
      }
      if (sumnorm.z() == 0.)
      {
        sumnorm += nZ;
      }
    }
  }
  if (noSurfaces == 0)
  {
#ifdef UDEBUG
    UUtils::Exception("USphere::SurfaceNormal(p)", "GeomSolids1002",
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
  return (noSurfaces > 0);
}


/////////////////////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

UVector3 USphere::ApproxSurfaceNormal(const UVector3& p) const
{
  ENorm side;
  UVector3 norm;
  double rho, rho2, radius, pPhi, pTheta;
  double distRMin, distRMax, distSPhi, distEPhi,
         distSTheta, distETheta, distMin;

  rho2 = p.x() * p.x() + p.y() * p.y();
  radius = std::sqrt(rho2 + p.z() * p.z());
  rho = std::sqrt(rho2);

  //
  // Distance to r shells
  //

  distRMax = std::fabs(radius - fRmax);
  if (fRmin)
  {
    distRMin = std::fabs(radius - fRmin);

    if (distRMin < distRMax)
    {
      distMin = distRMin;
      side = kNRMin;
    }
    else
    {
      distMin = distRMax;
      side = kNRMax;
    }
  }
  else
  {
    distMin = distRMax;
    side = kNRMax;
  }
  //
  // Distance to phi planes
  //
  // Protected against (0,0,z)

  pPhi = std::atan2(p.y(), p.x());
  if (pPhi < 0)
  {
    pPhi += 2 * UUtils::kPi;
  }

  if (!fFullPhiSphere && rho)
  {
    if (fSPhi < 0)
    {
      distSPhi = std::fabs(pPhi - (fSPhi + 2 * UUtils::kPi)) * rho;
    }
    else
    {
      distSPhi = std::fabs(pPhi - fSPhi) * rho;
    }

    distEPhi = std::fabs(pPhi - fSPhi - fDPhi) * rho;

    // Find new minimum
    //
    if (distSPhi < distEPhi)
    {
      if (distSPhi < distMin)
      {
        distMin = distSPhi;
        side = kNSPhi;
      }
    }
    else
    {
      if (distEPhi < distMin)
      {
        distMin = distEPhi;
        side = kNEPhi;
      }
    }
  }

  //
  // Distance to theta planes
  //

  if (!fFullThetaSphere && radius)
  {
    pTheta = std::atan2(rho, p.z());
    distSTheta = std::fabs(pTheta - fSTheta) * radius;
    distETheta = std::fabs(pTheta - fSTheta - fDTheta) * radius;

    // Find new minimum
    //
    if (distSTheta < distETheta)
    {
      if (distSTheta < distMin)
      {
        distMin = distSTheta;
        side = kNSTheta;
      }
    }
    else
    {
      if (distETheta < distMin)
      {
        distMin = distETheta;
        side = kNETheta;
      }
    }
  }

  switch (side)
  {
    case kNRMin:      // Inner radius
      norm = UVector3(-p.x() / radius, -p.y() / radius, -p.z() / radius);
      break;
    case kNRMax:      // Outer radius
      norm = UVector3(p.x() / radius, p.y() / radius, p.z() / radius);
      break;
    case kNSPhi:
      norm = UVector3(sinSPhi, -cosSPhi, 0);
      break;
    case kNEPhi:
      norm = UVector3(-sinEPhi, cosEPhi, 0);
      break;
    case kNSTheta:
      norm = UVector3(-cosSTheta * std::cos(pPhi),
                      -cosSTheta * std::sin(pPhi),
                      sinSTheta);
      break;
    case kNETheta:
      norm = UVector3(cosETheta * std::cos(pPhi),
                      cosETheta * std::sin(pPhi),
                      sinETheta);
      break;
    default:          // Should never reach this case ...

      UUtils::Exception("USphere::ApproxSurfaceNormal()",
                        "GeomSolids1002", UWarning, 1,
                        "Undefined side for valid surface normal to solid.");
      break;
  }

  return norm;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return UUtils::Infinity() if no intersection, or intersection distance <= tolerance
//
// -> If point is outside outer radius, compute intersection with rmax
//        - if no intersection return
//        - if  valid phi,theta return intersection Dist
//
// -> If shell, compute intersection with inner radius, taking largest +ve root
//        - if valid phi,theta, save intersection
//
// -> If phi segmented, compute intersection with phi half planes
//        - if valid intersection(r,theta), return smallest intersection of
//          inner shell & phi intersection
//
// -> If theta segmented, compute intersection with theta cones
//        - if valid intersection(r,phi), return smallest intersection of
//          inner shell & theta intersection
//
//
// NOTE:
// - `if valid' (above) implies tolerant checking of intersection points
//
// OPT:
// Move tolIO/ORmin/RMax2 precalcs to where they are needed -
// not required for most cases.
// Avoid atan2 for non theta cut USphere.

double USphere::DistanceToIn(const UVector3& p, const UVector3& v, double /*aPstep*/) const
{
  double snxt = UUtils::Infinity();     // snxt = default return value
  double rho2, rad2, pDotV2d, pDotV3d, pTheta;
  double tolSTheta = 0., tolETheta = 0.;
  const double dRmax = 100.*fRmax;

  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfAngTolerance = kAngTolerance * 0.5;
  const double halfTolerance = kTolerance * 0.5;
  const double halfRminTolerance = fRminTolerance * 0.5;
  const double tolORMin2 = (fRmin > halfRminTolerance)
                           ? (fRmin - halfRminTolerance) * (fRmin - halfRminTolerance) : 0;
  const double tolIRMin2 =
    (fRmin + halfRminTolerance) * (fRmin + halfRminTolerance);
  const double tolORMax2 =
    (fRmax + halfTolerance) * (fRmax + halfTolerance);
  const double tolIRMax2 =
    (fRmax - halfTolerance) * (fRmax - halfTolerance);

  // Intersection point
  //
  double xi, yi, zi, rhoi, rhoi2, radi2, iTheta;

  // Phi intersection
  //
  double Comp;

  // Phi precalcs
  //
  double Dist, cosPsi;

  // Theta precalcs
  //
  double dist2STheta, dist2ETheta;
  double t1, t2, b, c, d2, d, sd = UUtils::Infinity();

  // General Precalcs
  //
  rho2 = p.x() * p.x() + p.y() * p.y();
  rad2 = rho2 + p.z() * p.z();
  pTheta = std::atan2(std::sqrt(rho2), p.z());

  pDotV2d = p.x() * v.x() + p.y() * v.y();
  pDotV3d = pDotV2d + p.z() * v.z();

  // Theta precalcs
  //
  if (!fFullThetaSphere)
  {
    tolSTheta = fSTheta - halfAngTolerance;
    tolETheta = eTheta + halfAngTolerance;
  }

  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer USphere heading through solid (-> 0)
  // - No intersect -> no intersection with USphere
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2sd(pxvx+pyvy+pzvz)+sd^2(vx^2+vy^2+vz^2)=R^2
  // =>     rad2        +2sd(pDotV3d)      +sd^2                =R^2
  //
  // => sd=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  c = rad2 - fRmax * fRmax;

  if (c > kTolerance * fRmax)
  {
    // If outside tolerant boundary of outer USphere
    // [should be std::sqrt(rad2)-fRmax > halfTolerance]

    d2 = pDotV3d * pDotV3d - c;

    if (d2 >= 0)
    {
      sd = -pDotV3d - std::sqrt(d2);

      if (sd >= 0)
      {
        if (sd > dRmax) // Avoid rounding errors due to precision issues seen on
        {
          // 64 bits systems. Split long distances and recompute
          double fTerm = sd - std::fmod(sd, dRmax);
          sd = fTerm + DistanceToIn(p + fTerm * v, v);
        }
        xi   = p.x() + sd * v.x();
        yi   = p.y() + sd * v.y();
        rhoi = std::sqrt(xi * xi + yi * yi);

        if (!fFullPhiSphere && rhoi)    // Check phi intersection
        {
          cosPsi = (xi * cosCPhi + yi * sinCPhi) / rhoi;

          if (cosPsi >= cosHDPhiOT)
          {
            if (!fFullThetaSphere)   // Check theta intersection
            {
              zi = p.z() + sd * v.z();

              // rhoi & zi can never both be 0
              // (=>intersect at origin =>fRmax=0)
              //
              iTheta = std::atan2(rhoi, zi);
              if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
              {
                return snxt = sd;
              }
            }
            else
            {
              return snxt = sd;
            }
          }
        }
        else
        {
          if (!fFullThetaSphere)    // Check theta intersection
          {
            zi = p.z() + sd * v.z();

            // rhoi & zi can never both be 0
            // (=>intersect at origin => fRmax=0 !)
            //
            iTheta = std::atan2(rhoi, zi);
            if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
            {
              return snxt = sd;
            }
          }
          else
          {
            return snxt = sd;
          }
        }
      }
    }
    else    // No intersection with USphere
    {
      return snxt = UUtils::Infinity();
    }
  }
  else
  {
    // Inside outer radius
    // check not inside, and heading through USphere (-> 0 to in)

    d2 = pDotV3d * pDotV3d - c;

    if ((rad2 > tolIRMax2)
        && ((d2 >= kTolerance * fRmax) && (pDotV3d < 0)))
    {
      if (!fFullPhiSphere)
      {
        // Use inner phi tolerant boundary -> if on tolerant
        // phi boundaries, phi intersect code handles leaving/entering checks

        cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(rho2);

        if (cosPsi >= cosHDPhiIT)
        {
          // inside radii, delta r -ve, inside phi

          if (!fFullThetaSphere)
          {
            if ((pTheta >= tolSTheta + kAngTolerance)
                && (pTheta <= tolETheta - kAngTolerance))
            {
              return snxt = 0;
            }
          }
          else    // strictly inside Theta in both cases
          {
            return snxt = 0;
          }
        }
      }
      else
      {
        if (!fFullThetaSphere)
        {
          if ((pTheta >= tolSTheta + kAngTolerance)
              && (pTheta <= tolETheta - kAngTolerance))
          {
            return snxt = 0;
          }
        }
        else   // strictly inside Theta in both cases
        {
          return snxt = 0;
        }
      }
    }
  }

  // Inner spherical shell intersection
  // - Always farthest root, because would have passed through outer
  //   surface first.
  // - Tolerant check if travelling through solid

  if (fRmin)
  {
    c = rad2 - fRmin * fRmin;
    d2 = pDotV3d * pDotV3d - c;

    // Within tolerance inner radius of inner USphere
    // Check for immediate entry/already inside and travelling outwards

    if ((c > -halfRminTolerance) && (rad2 < tolIRMin2)
        && ((d2 < fRmin * VUSolid::Tolerance()) || (pDotV3d >= 0)))
    {
      if (!fFullPhiSphere)
      {
        // Use inner phi tolerant boundary -> if on tolerant
        // phi boundaries, phi intersect code handles leaving/entering checks

        cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(rho2);
        if (cosPsi >= cosHDPhiIT)
        {
          // inside radii, delta r -ve, inside phi
          //
          if (!fFullThetaSphere)
          {
            if ((pTheta >= tolSTheta + kAngTolerance)
                && (pTheta <= tolETheta - kAngTolerance))
            {
              return snxt = 0;
            }
          }
          else
          {
            return snxt = 0;
          }
        }
      }
      else
      {
        if (!fFullThetaSphere)
        {
          if ((pTheta >= tolSTheta + kAngTolerance)
              && (pTheta <= tolETheta - kAngTolerance))
          {
            return snxt = 0;
          }
        }
        else
        {
          return snxt = 0;
        }
      }
    }
    else   // Not special tolerant case
    {
      if (d2 >= 0)
      {
        sd = -pDotV3d + std::sqrt(d2);
        if (sd >= halfRminTolerance)  // It was >= 0 ??
        {
          xi   = p.x() + sd * v.x();
          yi   = p.y() + sd * v.y();
          rhoi = std::sqrt(xi * xi + yi * yi);

          if (!fFullPhiSphere && rhoi)   // Check phi intersection
          {
            cosPsi = (xi * cosCPhi + yi * sinCPhi) / rhoi;

            if (cosPsi >= cosHDPhiOT)
            {
              if (!fFullThetaSphere)  // Check theta intersection
              {
                zi = p.z() + sd * v.z();

                // rhoi & zi can never both be 0
                // (=>intersect at origin =>fRmax=0)
                //
                iTheta = std::atan2(rhoi, zi);
                if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
          else
          {
            if (!fFullThetaSphere)   // Check theta intersection
            {
              zi = p.z() + sd * v.z();

              // rhoi & zi can never both be 0
              // (=>intersect at origin => fRmax=0 !)
              //
              iTheta = std::atan2(rhoi, zi);
              if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
              {
                snxt = sd;
              }
            }
            else
            {
              snxt = sd;
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
  //
  if (!fFullPhiSphere)
  {
    // First phi surface ('S'tarting phi)
    // Comp = Component in outwards normal dirn
    //
    Comp = v.x() * sinSPhi - v.y() * cosSPhi;

    if (Comp < 0)
    {
      Dist = p.y() * cosSPhi - p.x() * sinSPhi;

      if (Dist < halfCarTolerance)
      {
        sd = Dist / Comp;

        if (sd < snxt)
        {
          if (sd > 0)
          {
            xi = p.x() + sd * v.x();
            yi = p.y() + sd * v.y();
            zi = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi ;
            radi2 = rhoi2 + zi * zi ;
          }
          else
          {
            sd    = 0;
            xi    = p.x();
            yi    = p.y();
            zi    = p.z();
            rhoi2 = rho2;
            radi2 = rad2;
          }
          if ((radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && ((yi * cosCPhi - xi * sinCPhi) <= 0))
          {
            // Check theta intersection
            // rhoi & zi can never both be 0
            // (=>intersect at origin =>fRmax=0)
            //
            if (!fFullThetaSphere)
            {
              iTheta = std::atan2(std::sqrt(rhoi2), zi);
              if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
              {
                // r and theta intersections good
                // - check intersecting with correct half-plane

                if ((yi * cosCPhi - xi * sinCPhi) <= 0)
                {
                  snxt = sd;
                }
              }
            }
            else
            {
              snxt = sd;
            }
          }
        }
      }
    }

    // Second phi surface ('E'nding phi)
    // Component in outwards normal dirn

    Comp = -(v.x() * sinEPhi - v.y() * cosEPhi);

    if (Comp < 0)
    {
      Dist = -(p.y() * cosEPhi - p.x() * sinEPhi);
      if (Dist < halfCarTolerance)
      {
        sd = Dist / Comp;

        if (sd < snxt)
        {
          if (sd > 0)
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi ;
            radi2 = rhoi2 + zi * zi ;
          }
          else
          {
            sd    = 0   ;
            xi    = p.x();
            yi    = p.y();
            zi    = p.z();
            rhoi2 = rho2  ;
            radi2 = rad2  ;
          }
          if ((radi2 <= tolORMax2)
              && (radi2 >= tolORMin2)
              && ((yi * cosCPhi - xi * sinCPhi) >= 0))
          {
            // Check theta intersection
            // rhoi & zi can never both be 0
            // (=>intersect at origin =>fRmax=0)
            //
            if (!fFullThetaSphere)
            {
              iTheta = std::atan2(std::sqrt(rhoi2), zi);
              if ((iTheta >= tolSTheta) && (iTheta <= tolETheta))
              {
                // r and theta intersections good
                // - check intersecting with correct half-plane

                if ((yi * cosCPhi - xi * sinCPhi) >= 0)
                {
                  snxt = sd;
                }
              }
            }
            else
            {
              snxt = sd;
            }
          }
        }
      }
    }
  }

  // Theta segment intersection

  if (!fFullThetaSphere)
  {

    // Intersection with theta surfaces
    // Known failure cases:
    // o  Inside tolerance of stheta surface, skim
    //    ~parallel to cone and Hit & enter etheta surface [& visa versa]
    //
    //    To solve: Check 2nd root of etheta surface in addition to stheta
    //
    // o  start/end theta is exactly UUtils::kPi/2
    // Intersections with cones
    //
    // Cone equation: x^2+y^2=z^2tan^2(t)
    //
    // => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
    //
    // => (px^2+py^2-pz^2tan^2(t))+2sd(pxvx+pyvy-pzvztan^2(t))
    //       + sd^2(vx^2+vy^2-vz^2tan^2(t)) = 0
    //
    // => sd^2(1-vz^2(1+tan^2(t))+2sd(pdotv2d-pzvztan^2(t))+(rho2-pz^2tan^2(t))=0

    if (fSTheta)
    {
      dist2STheta = rho2 - p.z() * p.z() * tanSTheta2;
    }
    else
    {
      dist2STheta = UUtils::Infinity();
    }
    if (eTheta < UUtils::kPi)
    {
      dist2ETheta = rho2 - p.z() * p.z() * tanETheta2;
    }
    else
    {
      dist2ETheta = UUtils::Infinity();
    }
    if (pTheta < tolSTheta)
    {
      // Inside (theta<stheta-tol) stheta cone
      // First root of stheta cone, second if first root -ve

      t1 = 1 - v.z() * v.z() * (1 + tanSTheta2);
      t2 = pDotV2d - p.z() * v.z() * tanSTheta2;
      if (t1)
      {
        b = t2 / t1;
        c = dist2STheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b - d;    // First root
          zi = p.z() + sd * v.z();

          if ((sd < 0) || (zi * (fSTheta - UUtils::kPi / 2) > 0))
          {
            sd = -b + d;  // Second root
          }
          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;
            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (fSTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2) // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }

      // Possible intersection with ETheta cone.
      // Second >= 0 root should be considered

      if (eTheta < UUtils::kPi)
      {
        t1 = 1 - v.z() * v.z() * (1 + tanETheta2);
        t2 = pDotV2d - p.z() * v.z() * tanETheta2;
        if (t1)
        {
          b = t2 / t1;
          c = dist2ETheta / t1;
          d2 = b * b - c;

          if (d2 >= 0)
          {
            d = std::sqrt(d2);
            sd = -b + d;    // Second root

            if ((sd >= 0) && (sd < snxt))
            {
              xi    = p.x() + sd * v.x();
              yi    = p.y() + sd * v.y();
              zi    = p.z() + sd * v.z();
              rhoi2 = xi * xi + yi * yi;
              radi2 = rhoi2 + zi * zi;

              if ((radi2 <= tolORMax2)
                  && (radi2 >= tolORMin2)
                  && (zi * (eTheta - UUtils::kPi / 2) <= 0))
              {
                if (!fFullPhiSphere && rhoi2)  // Check phi intersection
                {
                  cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                  if (cosPsi >= cosHDPhiOT)
                  {
                    snxt = sd;
                  }
                }
                else
                {
                  snxt = sd;
                }
              }
            }
          }
        }
      }
    }
    else if (pTheta > tolETheta)
    {
      // dist2ETheta<-kRadTolerance*0.5 && dist2STheta>0)
      // Inside (theta > etheta+tol) e-theta cone
      // First root of etheta cone, second if first root 'imaginary'

      t1 = 1 - v.z() * v.z() * (1 + tanETheta2);
      t2 = pDotV2d - p.z() * v.z() * tanETheta2;
      if (t1)
      {
        b = t2 / t1;
        c = dist2ETheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b - d;    // First root
          zi = p.z() + sd * v.z();

          if ((sd < 0) || (zi * (eTheta - UUtils::kPi / 2) > 0))
          {
            sd = -b + d;           // second root
          }
          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;

            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (eTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2) // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }

      // Possible intersection with STheta cone.
      // Second >= 0 root should be considered

      if (fSTheta)
      {
        t1 = 1 - v.z() * v.z() * (1 + tanSTheta2);
        t2 = pDotV2d - p.z() * v.z() * tanSTheta2;
        if (t1)
        {
          b = t2 / t1;
          c = dist2STheta / t1;
          d2 = b * b - c;

          if (d2 >= 0)
          {
            d = std::sqrt(d2);
            sd = -b + d;    // Second root

            if ((sd >= 0) && (sd < snxt))
            {
              xi    = p.x() + sd * v.x();
              yi    = p.y() + sd * v.y();
              zi    = p.z() + sd * v.z();
              rhoi2 = xi * xi + yi * yi;
              radi2 = rhoi2 + zi * zi;

              if ((radi2 <= tolORMax2)
                  && (radi2 >= tolORMin2)
                  && (zi * (fSTheta - UUtils::kPi / 2) <= 0))
              {
                if (!fFullPhiSphere && rhoi2)  // Check phi intersection
                {
                  cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                  if (cosPsi >= cosHDPhiOT)
                  {
                    snxt = sd;
                  }
                }
                else
                {
                  snxt = sd;
                }
              }
            }
          }
        }
      }
    }
    else if ((pTheta < tolSTheta + kAngTolerance)
             && (fSTheta > halfAngTolerance))
    {
      // In tolerance of stheta
      // If entering through solid [r,phi] => 0 to in
      // else try 2nd root

      t2 = pDotV2d - p.z() * v.z() * tanSTheta2;
      if ((t2 >= 0 && tolIRMin2 < rad2 && rad2 < tolIRMax2 && fSTheta < UUtils::kPi / 2)
          || (t2 < 0  && tolIRMin2 < rad2 && rad2 < tolIRMax2 && fSTheta > UUtils::kPi / 2)
          || (v.z() < 0 && tolIRMin2 < rad2 && rad2 < tolIRMax2 && fSTheta == UUtils::kPi / 2))
      {
        if (!fFullPhiSphere && rho2)  // Check phi intersection
        {
          cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(rho2);
          if (cosPsi >= cosHDPhiIT)
          {
            return 0;
          }
        }
        else
        {
          return 0;
        }
      }

      // Not entering immediately/travelling through

      t1 = 1 - v.z() * v.z() * (1 + tanSTheta2);
      if (t1)
      {
        b = t2 / t1;
        c = dist2STheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b + d;
          if ((sd >= halfCarTolerance) && (sd < snxt) && (fSTheta < UUtils::kPi / 2))
          {
            // ^^^^^^^^^^^^^^^^^^^^^ shouldn't it be >=0 instead ?
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;

            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (fSTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2)   // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }
    else if ((pTheta > tolETheta - kAngTolerance) && (eTheta < UUtils::kPi - kAngTolerance))
    {

      // In tolerance of etheta
      // If entering through solid [r,phi] => 0 to in
      // else try 2nd root

      t2 = pDotV2d - p.z() * v.z() * tanETheta2;

      if (((t2 < 0) && (eTheta < UUtils::kPi / 2)
           && (tolIRMin2 < rad2) && (rad2 < tolIRMax2))
          || ((t2 >= 0) && (eTheta > UUtils::kPi / 2)
              && (tolIRMin2 < rad2) && (rad2 < tolIRMax2))
          || ((v.z() > 0) && (eTheta == UUtils::kPi / 2)
              && (tolIRMin2 < rad2) && (rad2 < tolIRMax2)))
      {
        if (!fFullPhiSphere && rho2)   // Check phi intersection
        {
          cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / std::sqrt(rho2);
          if (cosPsi >= cosHDPhiIT)
          {
            return 0;
          }
        }
        else
        {
          return 0;
        }
      }

      // Not entering immediately/travelling through

      t1 = 1 - v.z() * v.z() * (1 + tanETheta2);
      if (t1)
      {
        b = t2 / t1;
        c = dist2ETheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b + d;

          if ((sd >= halfCarTolerance)
              && (sd < snxt) && (eTheta > UUtils::kPi / 2))
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;

            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (eTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2)  // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }
    else
    {
      // stheta+tol<theta<etheta-tol
      // For BOTH stheta & etheta check 2nd root for validity [r,phi]

      t1 = 1 - v.z() * v.z() * (1 + tanSTheta2);
      t2 = pDotV2d - p.z() * v.z() * tanSTheta2;
      if (t1)
      {
        b = t2 / t1;
        c = dist2STheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b + d;    // second root

          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;

            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (fSTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2)  // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
      t1 = 1 - v.z() * v.z() * (1 + tanETheta2);
      t2 = pDotV2d - p.z() * v.z() * tanETheta2;
      if (t1)
      {
        b = t2 / t1;
        c = dist2ETheta / t1;
        d2 = b * b - c;

        if (d2 >= 0)
        {
          d = std::sqrt(d2);
          sd = -b + d;    // second root

          if ((sd >= 0) && (sd < snxt))
          {
            xi    = p.x() + sd * v.x();
            yi    = p.y() + sd * v.y();
            zi    = p.z() + sd * v.z();
            rhoi2 = xi * xi + yi * yi;
            radi2 = rhoi2 + zi * zi;

            if ((radi2 <= tolORMax2)
                && (radi2 >= tolORMin2)
                && (zi * (eTheta - UUtils::kPi / 2) <= 0))
            {
              if (!fFullPhiSphere && rhoi2)  // Check phi intersection
              {
                cosPsi = (xi * cosCPhi + yi * sinCPhi) / std::sqrt(rhoi2);
                if (cosPsi >= cosHDPhiOT)
                {
                  snxt = sd;
                }
              }
              else
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to radial planes
// - Only to phi planes if outside phi extent
// - Only to theta planes if outside theta extent
// - Return 0 if point inside

double USphere::SafetyFromOutside(const UVector3& p, bool /*aAccurate*/) const
{
  double safe = 0.0, safeRMin, safeRMax, safePhi, safeTheta;
  double rho2, rds, rho;
  double cosPsi;
  double pTheta, dTheta1, dTheta2;
  rho2 = p.x() * p.x() + p.y() * p.y();
  rds = std::sqrt(rho2 + p.z() * p.z());
  rho = std::sqrt(rho2);

  //
  // Distance to r shells
  //
  if (fRmin)
  {
    safeRMin = fRmin - rds;
    safeRMax = rds - fRmax;
    if (safeRMin > safeRMax)
    {
      safe = safeRMin;
    }
    else
    {
      safe = safeRMax;
    }
  }
  else
  {
    safe = rds - fRmax;
  }

  //
  // Distance to phi extent
  //
  if (!fFullPhiSphere && rho)
  {
    // Psi=angle from central phi to point
    //
    cosPsi = (p.x() * cosCPhi + p.y() * sinCPhi) / rho;
    if (cosPsi < std::cos(hDPhi))
    {
      // Point lies outside phi range
      //
      if ((p.y() * cosCPhi - p.x() * sinCPhi) <= 0)
      {
        safePhi = std::fabs(p.x() * sinSPhi - p.y() * cosSPhi);
      }
      else
      {
        safePhi = std::fabs(p.x() * sinEPhi - p.y() * cosEPhi);
      }
      if (safePhi > safe)
      {
        safe = safePhi;
      }
    }
  }
  //
  // Distance to Theta extent
  //
  if ((rds != 0.0) && (!fFullThetaSphere))
  {
    pTheta = std::acos(p.z() / rds);
    if (pTheta < 0)
    {
      pTheta += UUtils::kPi;
    }
    dTheta1 = fSTheta - pTheta;
    dTheta2 = pTheta - eTheta;
    if (dTheta1 > dTheta2)
    {
      if (dTheta1 >= 0)          // WHY ???????????
      {
        safeTheta = rds * std::sin(dTheta1);
        if (safe <= safeTheta)
        {
          safe = safeTheta;
        }
      }
    }
    else
    {
      if (dTheta2 >= 0)
      {
        safeTheta = rds * std::sin(dTheta2);
        if (safe <= safeTheta)
        {
          safe = safeTheta;
        }
      }
    }
  }

  if (safe < 0)
  {
    safe = 0;
  }
  return safe;
}

/////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from 'inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

double USphere::DistanceToOut(const UVector3& p, const UVector3& v, UVector3& n, bool& validNorm, double /*aPstep*/) const
{
  double snxt = UUtils::Infinity();    // snxt is default return value
  double sphi = UUtils::Infinity(), stheta = UUtils::Infinity();
  ESide side = kNull, sidephi = kNull, sidetheta = kNull;

  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfAngTolerance = kAngTolerance * 0.5;
  const double halfTolerance = kTolerance * 0.5;
  const double halfRminTolerance = fRminTolerance * 0.5;
  const double rMaxPlus = fRmax + halfTolerance;
  const double Rmin_minus = (fRmin) ? fRmin - halfRminTolerance : 0;
  double t1, t2;
  double b, c, d;

  // Variables for phi intersection:

  double pDistS, compS, pDistE, compE, sphi2, vphi;

  double rho2, rad2, pDotV2d, pDotV3d;

  double xi, yi, zi;    // Intersection point

  // Theta precals
  //
  double rhoSecTheta;
  double dist2STheta, dist2ETheta, distTheta;
  double d2, sd;

  // General Precalcs
  //
  rho2 = p.x() * p.x() + p.y() * p.y();
  rad2 = rho2 + p.z() * p.z();

  pDotV2d = p.x() * v.x() + p.y() * v.y();
  pDotV3d = pDotV2d + p.z() * v.z();

  // Radial Intersections from USphere::DistanceToIn
  //
  // Outer spherical shell intersection
  // - Only if outside tolerant fRmax
  // - Check for if inside and outer USphere heading through solid (-> 0)
  // - No intersect -> no intersection with USphere
  //
  // Shell eqn: x^2+y^2+z^2=RSPH^2
  //
  // => (px+svx)^2+(py+svy)^2+(pz+svz)^2=R^2
  //
  // => (px^2+py^2+pz^2) +2sd(pxvx+pyvy+pzvz)+sd^2(vx^2+vy^2+vz^2)=R^2
  // =>     rad2        +2sd(pDotV3d)      +sd^2                =R^2
  //
  // => sd=-pDotV3d+-std::sqrt(pDotV3d^2-(rad2-R^2))

  if ((rad2 <= rMaxPlus * rMaxPlus) && (rad2 >= Rmin_minus * Rmin_minus))
  {
    c = rad2 - fRmax * fRmax;

    if (c < kTolerance * fRmax)
    {
      // Within tolerant Outer radius
      //
      // The test is
      //     rad  - fRmax < 0.5*kRadTolerance
      // => rad < fRmax + 0.5*kRadTol
      // => rad2 < (fRmax + 0.5*kRadTol)^2
      // => rad2 < fRmax^2 + 2.*0.5*fRmax*kRadTol + 0.25*kRadTol*kRadTol
      // => rad2 - fRmax^2    <~    fRmax*kRadTol

      d2 = pDotV3d * pDotV3d - c;

      if ((c > - kTolerance * fRmax)   // on tolerant surface
          && ((pDotV3d >= 0) || (d2 < 0)))    // leaving outside from Rmax
        // not re-entering
      {
        validNorm = true;
        n        = UVector3(p.x() / fRmax, p.y() / fRmax, p.z() / fRmax);
        return snxt = 0;
      }
      else
      {
        snxt = -pDotV3d + std::sqrt(d2);  // second root since inside Rmax
        side =  kRMax;
      }
    }

    // Inner spherical shell intersection:
    // Always first >=0 root, because would have passed
    // from outside of Rmin surface .

    if (fRmin)
    {
      c = rad2 - fRmin * fRmin;
      d2 = pDotV3d * pDotV3d - c;

      if (c > - fRminTolerance * fRmin) // 2.0 * (0.5*kRadTolerance) * fRmin
      {
        if ((c < fRminTolerance * fRmin)            // leaving from Rmin
            && (d2 >= fRminTolerance * fRmin) && (pDotV3d < 0))
        {
          validNorm = false;  // Rmin surface is concave
          n        = UVector3(-p.x() / fRmin, -p.y() / fRmin, -p.z() / fRmin);
          return snxt = 0;
        }
        else
        {
          if (d2 >= 0.)
          {
            sd = -pDotV3d - std::sqrt(d2);

            if (sd >= 0.)    // Always intersect Rmin first
            {
              snxt = sd;
              side = kRMin;
            }
          }
        }
      }
    }
  }

  // Theta segment intersection

  if (!fFullThetaSphere)
  {
    // Intersection with theta surfaces
    //
    // Known failure cases:
    // o  Inside tolerance of stheta surface, skim
    //    ~parallel to cone and Hit & enter etheta surface [& visa versa]
    //
    //    To solve: Check 2nd root of etheta surface in addition to stheta
    //
    // o  start/end theta is exactly UUtils::kPi/2
    //
    // Intersections with cones
    //
    // Cone equation: x^2+y^2=z^2tan^2(t)
    //
    // => (px+svx)^2+(py+svy)^2=(pz+svz)^2tan^2(t)
    //
    // => (px^2+py^2-pz^2tan^2(t))+2sd(pxvx+pyvy-pzvztan^2(t))
    //       + sd^2(vx^2+vy^2-vz^2tan^2(t)) = 0
    //
    // => sd^2(1-vz^2(1+tan^2(t))+2sd(pdotv2d-pzvztan^2(t))+(rho2-pz^2tan^2(t))=0
    //

    if (fSTheta) // intersection with first cons
    {
      if (std::fabs(tanSTheta) > 5. / kAngTolerance) // kons is plane z=0
      {
        if (v.z() > 0.)
        {
          if (std::fabs(p.z()) <= halfTolerance)
          {
            validNorm = true;
            n = UVector3(0., 0., 1.);
            return snxt = 0;
          }
          stheta    = -p.z() / v.z();
          sidetheta = kSTheta;
        }
      }
      else // kons is not plane
      {
        t1          = 1 - v.z() * v.z() * (1 + tanSTheta2);
        t2          = pDotV2d - p.z() * v.z() * tanSTheta2; // ~vDotN if p on cons
        dist2STheta = rho2 - p.z() * p.z() * tanSTheta2; // t3

        distTheta = std::sqrt(rho2) - p.z() * tanSTheta;

        if (std::fabs(t1) < halfAngTolerance) // 1st order equation,
        {
          // v parallel to kons
          if (v.z() > 0.)
          {
            if (std::fabs(distTheta) < halfTolerance) // p on surface
            {
              if ((fSTheta < UUtils::kPi / 2) && (p.z() > 0.))
              {
                validNorm = false;
                if (rho2)
              {
                rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));

                n = UVector3(p.x() / rhoSecTheta,
                             p.y() / rhoSecTheta,
                             -std::sin(fSTheta));
              }
              else
              {
                n = UVector3(0., 0., 1.);
              }
                return snxt = 0.;
              }
              else if ((fSTheta > UUtils::kPi / 2) && (p.z() <= 0))
              {
                validNorm = true;
                if (rho2)
                {
                  rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));

                  n = UVector3(p.x() / rhoSecTheta,
                               p.y() / rhoSecTheta,
                               std::sin(fSTheta));
                }
                else n = UVector3(0., 0., 1.);
                return snxt = 0.;
              }
            }
            stheta    = -0.5 * dist2STheta / t2;
            sidetheta = kSTheta;
          }
        }     // 2nd order equation, 1st root of fSTheta cone,
        else   // 2nd if 1st root -ve
        {
          if (std::fabs(distTheta) < halfTolerance)
          {
            if ((fSTheta > UUtils::kPi / 2) && (t2 >= 0.)) // leave
            {
              validNorm = true;
              if (rho2)
              {
                rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));

                n = UVector3(p.x() / rhoSecTheta,
                             p.y() / rhoSecTheta,
                             std::sin(fSTheta));
              }
              else
              {
                n = UVector3(0., 0., 1.);
              }
              return snxt = 0.;
            }
            else if ((fSTheta < UUtils::kPi / 2) && (t2 < 0.) && (p.z() >= 0.)) // leave
            {
              validNorm = false;
              if (rho2)
              {
                rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));

                n = UVector3(p.x() / rhoSecTheta,
                             p.y() / rhoSecTheta,
                             -std::sin(fSTheta));
              }
              else
              {
                n = UVector3(0., 0., 1.);
              }

              return snxt = 0.;
            }
          }
          b = t2 / t1;
          c = dist2STheta / t1;
          d2 = b * b - c;

          if (d2 >= 0.)
          {
            d = std::sqrt(d2);

            if (fSTheta > UUtils::kPi / 2)
            {
              sd = -b - d;         // First root

              if (((std::fabs(sd) < halfTolerance) && (t2 < 0.))
                  || (sd < 0.) || ((sd > 0.) && (p.z() + sd * v.z() > 0.)))
              {
                sd = -b + d; // 2nd root
              }
              if ((sd > halfTolerance) && (p.z() + sd * v.z() <= 0.))
              {
                stheta    = sd;
                sidetheta = kSTheta;
              }
            }
            else // sTheta < UUtils::kPi/2, concave surface, no normal
            {
              sd = -b - d;         // First root

              if (((std::fabs(sd) < halfTolerance) && (t2 >= 0.))
                  || (sd < 0.) || ((sd > 0.) && (p.z() + sd * v.z() < 0.)))
              {
                sd = -b + d; // 2nd root
              }
              if ((sd > halfTolerance) && (p.z() + sd * v.z() >= 0.))
              {
                stheta    = sd;
                sidetheta = kSTheta;
              }
            }
          }
        }
      }
    }
    if (eTheta < UUtils::kPi) // intersection with second cons
    {
      if (std::fabs(tanETheta) > 5. / kAngTolerance) // kons is plane z=0
      {
        if (v.z() < 0.)
        {
          if (std::fabs(p.z()) <= halfTolerance)
          {
            validNorm = true;
            n = UVector3(0., 0., -1.);
            return snxt = 0;
          }
          sd = -p.z() / v.z();

          if (sd < stheta)
          {
            stheta    = sd;
            sidetheta = kETheta;
          }
        }
      }
      else // kons is not plane
      {
        t1          = 1 - v.z() * v.z() * (1 + tanETheta2);
        t2          = pDotV2d - p.z() * v.z() * tanETheta2; // ~vDotN if p on cons
        dist2ETheta = rho2 - p.z() * p.z() * tanETheta2; // t3

        distTheta = std::sqrt(rho2) - p.z() * tanETheta;

        if (std::fabs(t1) < halfAngTolerance) // 1st order equation,
        {
          // v parallel to kons
          if (v.z() < 0.)
          {
            if (std::fabs(distTheta) < halfTolerance) // p on surface
            {
              if ((eTheta > UUtils::kPi / 2) && (p.z() < 0.))
              {
                validNorm = false;
                 if (rho2)
                {
                  rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
                  n = UVector3(p.x() / rhoSecTheta,
                               p.y() / rhoSecTheta,
                               sinETheta);
                }
                else
                {
                  n = UVector3(0., 0., -1.);
                }
                return snxt = 0.;
              }
              else if ((eTheta < UUtils::kPi / 2) && (p.z() >= 0))
              {
                validNorm = true;
                if (rho2)
                {
                  rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
                  n = UVector3(p.x() / rhoSecTheta,
                               p.y() / rhoSecTheta,
                               -sinETheta);
                }
                else
                {
                  n = UVector3(0., 0., -1.);
                }
                return snxt = 0.;
              }
            }
            sd = -0.5 * dist2ETheta / t2;

            if (sd < stheta)
            {
              stheta    = sd;
              sidetheta = kETheta;
            }
          }
        }     // 2nd order equation, 1st root of fSTheta cone
        else   // 2nd if 1st root -ve
        {
          if (std::fabs(distTheta) < halfTolerance)
          {
            if ((eTheta < UUtils::kPi / 2) && (t2 >= 0.)) // leave
            {
              validNorm = true;
              if (rho2)
              {
                rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
                n = UVector3(p.x() / rhoSecTheta,
                             p.y() / rhoSecTheta,
                             -sinETheta);
              }
              else n = UVector3(0., 0., -1.);
              return snxt = 0.;
            }
            else if ((eTheta > UUtils::kPi / 2)
                     && (t2 < 0.) && (p.z() <= 0.)) // leave
            {
              validNorm = false;
               if (rho2)
              {
                rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
                n = -UVector3(p.x() / rhoSecTheta,
                             p.y() / rhoSecTheta,
                             sinETheta);
              }
              else n = UVector3(0., 0., -1.);
             
              return snxt = 0.;
            }
          }
          b = t2 / t1;
          c = dist2ETheta / t1;
          d2 = b * b - c;

          if (d2 >= 0.)
          {
            d = std::sqrt(d2);

            if (eTheta < UUtils::kPi / 2)
            {
              sd = -b - d;         // First root

              if (((std::fabs(sd) < halfTolerance) && (t2 < 0.))
                  || (sd < 0.))
              {
                sd = -b + d; // 2nd root
              }
              if (sd > halfTolerance)
              {
                if (sd < stheta)
                {
                  stheta    = sd;
                  sidetheta = kETheta;
                }
              }
            }
            else // sTheta+fDTheta > UUtils::kPi/2, concave surface, no normal
            {
              sd = -b - d;         // First root

              if (((std::fabs(sd) < halfTolerance) && (t2 >= 0.))
                  || (sd < 0.) || ((sd > 0.) && (p.z() + sd * v.z() > 0.)))
              {
                sd = -b + d; // 2nd root
              }
              if ((sd > halfTolerance) && (p.z() + sd * v.z() <= 0.))
              {
                if (sd < stheta)
                {
                  stheta    = sd;
                  sidetheta = kETheta;
                }
              }
            }
          }
        }
      }
    }

  } // end theta intersections

  // Phi Intersection

  if (!fFullPhiSphere)
  {
    if (p.x() || p.y()) // Check if on z axis (rho not needed later)
    {
      // pDist -ve when inside

      pDistS = p.x() * sinSPhi - p.y() * cosSPhi;
      pDistE = -p.x() * sinEPhi + p.y() * cosEPhi;

      // Comp -ve when in direction of outwards normal

      compS  = -sinSPhi * v.x() + cosSPhi * v.y();
      compE  =  sinEPhi * v.x() - cosEPhi * v.y();
      sidephi = kNull;

      if ((pDistS <= 0) && (pDistE <= 0))
      {
        // Inside both phi *full* planes

        if (compS < 0)
        {
          sphi = pDistS / compS;
          xi   = p.x() + sphi * v.x();
          yi   = p.y() + sphi * v.y();

          // Check intersection with correct half-plane (if not -> no intersect)
          //
          if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
          {
            vphi = std::atan2(v.y(), v.x());
            sidephi = kSPhi;
            if (((fSPhi - halfAngTolerance) <= vphi)
                && ((ePhi + halfAngTolerance) >= vphi))
            {
              sphi = UUtils::Infinity();
            }
          }
          else if ((yi * cosCPhi - xi * sinCPhi) >= 0)
          {
            sphi = UUtils::Infinity();
          }
          else
          {
            sidephi = kSPhi;
            if (pDistS > -halfCarTolerance)
            {
              sphi = 0;  // Leave by sphi
            }
          }
        }
        else
        {
          sphi = UUtils::Infinity();
        }

        if (compE < 0)
        {
          sphi2 = pDistE / compE;
          if (sphi2 < sphi) // Only check further if < starting phi intersection
          {
            xi = p.x() + sphi2 * v.x();
            yi = p.y() + sphi2 * v.y();

            // Check intersection with correct half-plane
            //
            if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
            {
              // Leaving via ending phi
              //
              vphi = std::atan2(v.y(), v.x());

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
            else if ((yi * cosCPhi - xi * sinCPhi) >= 0) // Leaving via ending phi
            {
              sidephi = kEPhi;
              if (pDistE <= -halfCarTolerance)
              {
                sphi = sphi2;
              }
              else
              {
                sphi = 0;
              }
            }
          }
        }
      }
      else if ((pDistS >= 0) && (pDistE >= 0)) // Outside both *full* phi planes
      {
        if (pDistS <= pDistE)
        {
          sidephi = kSPhi;
        }
        else
        {
          sidephi = kEPhi;
        }
        if (fDPhi > UUtils::kPi)
        {
          if ((compS < 0) && (compE < 0))
          {
            sphi = 0;
          }
          else
          {
            sphi = UUtils::Infinity();
          }
        }
        else
        {
          // if towards both >=0 then once inside (after error)
          // will remain inside

          if ((compS >= 0) && (compE >= 0))
          {
            sphi = UUtils::Infinity();
          }
          else
          {
            sphi = 0;
          }
        }
      }
      else if ((pDistS > 0) && (pDistE < 0))
      {
        // Outside full starting plane, inside full ending plane

        if (fDPhi > UUtils::kPi)
        {
          if (compE < 0)
          {
            sphi = pDistE / compE;
            xi   = p.x() + sphi * v.x();
            yi   = p.y() + sphi * v.y();

            // Check intersection in correct half-plane
            // (if not -> not leaving phi extent)
            //
            if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
            {
              vphi = std::atan2(v.y(), v.x());
              sidephi = kSPhi;
              if (((fSPhi - halfAngTolerance) <= vphi)
                  && ((ePhi + halfAngTolerance) >= vphi))
              {
                sphi = UUtils::Infinity();
              }
            }
            else if ((yi * cosCPhi - xi * sinCPhi) <= 0)
            {
              sphi = UUtils::Infinity();
            }
            else // Leaving via Ending phi
            {
              sidephi = kEPhi;
              if (pDistE > -halfCarTolerance)
              {
                sphi = 0.;
              }
            }
          }
          else
          {
            sphi = UUtils::Infinity();
          }
        }
        else
        {
          if (compS >= 0)
          {
            if (compE < 0)
            {
              sphi = pDistE / compE;
              xi   = p.x() + sphi * v.x();
              yi   = p.y() + sphi * v.y();

              // Check intersection in correct half-plane
              // (if not -> remain in extent)
              //
              if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
              {
                vphi = std::atan2(v.y(), v.x());
                sidephi = kSPhi;
                if (((fSPhi - halfAngTolerance) <= vphi)
                    && ((ePhi + halfAngTolerance) >= vphi))
                {
                  sphi = UUtils::Infinity();
                }
              }
              else if ((yi * cosCPhi - xi * sinCPhi) <= 0)
              {
                sphi = UUtils::Infinity();
              }
              else // otherwise leaving via Ending phi
              {
                sidephi = kEPhi;
              }
            }
            else sphi = UUtils::Infinity();
          }
          else // leaving immediately by starting phi
          {
            sidephi = kSPhi;
            sphi    = 0;
          }
        }
      }
      else
      {
        // Must be pDistS < 0 && pDistE > 0
        // Inside full starting plane, outside full ending plane

        if (fDPhi > UUtils::kPi)
        {
          if (compS < 0)
          {
            sphi = pDistS / compS;
            xi = p.x() + sphi * v.x();
            yi = p.y() + sphi * v.y();

            // Check intersection in correct half-plane
            // (if not -> not leaving phi extent)
            //
            if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
            {
              vphi = std::atan2(v.y(), v.x());
              sidephi = kSPhi;
              if (((fSPhi - halfAngTolerance) <= vphi)
                  && ((ePhi + halfAngTolerance) >= vphi))
              {
                sphi = UUtils::Infinity();
              }
            }
            else if ((yi * cosCPhi - xi * sinCPhi) >= 0)
            {
              sphi = UUtils::Infinity();
            }
            else  // Leaving via Starting phi
            {
              sidephi = kSPhi;
              if (pDistS > -halfCarTolerance)
              {
                sphi = 0;
              }
            }
          }
          else
          {
            sphi = UUtils::Infinity();
          }
        }
        else
        {
          if (compE >= 0)
          {
            if (compS < 0)
            {
              sphi = pDistS / compS;
              xi   = p.x() + sphi * v.x();
              yi   = p.y() + sphi * v.y();

              // Check intersection in correct half-plane
              // (if not -> remain in extent)
              //
              if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
              {
                vphi = std::atan2(v.y(), v.x());
                sidephi = kSPhi;
                if (((fSPhi - halfAngTolerance) <= vphi)
                    && ((ePhi + halfAngTolerance) >= vphi))
                {
                  sphi = UUtils::Infinity();
                }
              }
              else if ((yi * cosCPhi - xi * sinCPhi) >= 0)
              {
                sphi = UUtils::Infinity();
              }
              else // otherwise leaving via Starting phi
              {
                sidephi = kSPhi;
              }
            }
            else
            {
              sphi = UUtils::Infinity();
            }
          }
          else // leaving immediately by ending
          {
            sidephi = kEPhi;
            sphi    = 0   ;
          }
        }
      }
    }
    else
    {
      // On z axis + travel not || to z axis -> if phi of vector direction
      // within phi of shape, Step limited by rmax, else Step =0

      if (v.x() || v.y())
      {
        vphi = std::atan2(v.y(), v.x());
        if ((fSPhi - halfAngTolerance < vphi) && (vphi < ePhi + halfAngTolerance))
        {
          sphi = UUtils::Infinity();
        }
        else
        {
          sidephi = kSPhi; // arbitrary
          sphi    = 0;
        }
      }
      else  // travel along z - no phi intersection
      {
        sphi = UUtils::Infinity();
      }
    }
    if (sphi < snxt)  // Order intersecttions
    {
      snxt = sphi;
      side = sidephi;
    }
  }
  if (stheta < snxt) // Order intersections
  {
    snxt = stheta;
    side = sidetheta;
  }

  switch (side)
  {
    case kRMax:
      xi = p.x() + snxt * v.x();
      yi = p.y() + snxt * v.y();
      zi = p.z() + snxt * v.z();
      n = UVector3(xi / fRmax, yi / fRmax, zi / fRmax);
      validNorm = true;
      break;

    case kRMin:
      xi = p.x() + snxt * v.x();
      yi = p.y() + snxt * v.y();
      zi = p.z() + snxt * v.z();
      n = UVector3(-xi / fRmin, -yi / fRmin, -zi / fRmin);
      validNorm = false; // Rmin is concave
      break;

    case kSPhi:
      if (fDPhi <= UUtils::kPi)    // Normal to Phi-
      {
        n = UVector3(sinSPhi, -cosSPhi, 0);
        validNorm = true;
      }
      else
      {
        n = UVector3(sinSPhi, -cosSPhi, 0);
        validNorm = false;
      }
      break;

    case kEPhi:
      if (fDPhi <= UUtils::kPi)     // Normal to Phi+
      {
        n = UVector3(-sinEPhi, cosEPhi, 0);
        validNorm = true;
      }
      else
      {
        n = UVector3(-sinEPhi, cosEPhi, 0);
        validNorm = false;
      }
      break;

    case kSTheta:
      if (fSTheta == UUtils::kPi / 2)
      {
        n = UVector3(0., 0., 1.);
        validNorm = true;
      }
      else if (fSTheta > UUtils::kPi / 2)
      {
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        rho2 = xi * xi + yi * yi;
        if (rho2)
        {
          rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));
          n = UVector3(xi / rhoSecTheta, yi / rhoSecTheta,
                       -tanSTheta / std::sqrt(1 + tanSTheta2));
        }
        else
        {
          n = UVector3(0., 0., 1.);
        }
        validNorm = true;
      }
      else
      {
       xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        rho2 = xi * xi + yi * yi;
        if (rho2)
        {
          rhoSecTheta = std::sqrt(rho2 * (1 + tanSTheta2));
          n = UVector3(xi / rhoSecTheta, yi / rhoSecTheta,
                       -tanSTheta / std::sqrt(1 + tanSTheta2));
        }
        else
        {
          n = UVector3(0., 0., 1.);
        }
      
        validNorm = false;  // Concave STheta cone
      }
      break;

    case kETheta:
      if (eTheta == UUtils::kPi / 2)
      {
        n        = UVector3(0., 0., -1.);
        validNorm = true;
      }
      else if (eTheta < UUtils::kPi / 2)
      {
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        rho2 = xi * xi + yi * yi;
        if (rho2)
        {
          rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
          n = UVector3(xi / rhoSecTheta, yi / rhoSecTheta,
                       -tanETheta / std::sqrt(1 + tanETheta2));
        }
        else
        {
          n = UVector3(0., 0., -1.);
        }
        validNorm = true;
      }
      else
      {
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        rho2 = xi * xi + yi * yi;
        if (rho2)
        {
          rhoSecTheta = std::sqrt(rho2 * (1 + tanETheta2));
          n = -UVector3(xi / rhoSecTheta, yi / rhoSecTheta,
                       -tanETheta / std::sqrt(1 + tanETheta2));
        }
        else
        {
          n = UVector3(0., 0., -1.);
        }
        validNorm = false;  // Concave ETheta cone
      }
      break;

    default:
      cout << std::endl;

      std::ostringstream message;
      int oldprc = message.precision(16);
      message << "Undefined side for valid surface normal to solid."
              << std::endl
              << "Position:"  << std::endl << std::endl
              << "p.x = "  << p.x() << " mm" << std::endl
              << "p.y = "  << p.y() << " mm" << std::endl
              << "p.z = "  << p.z() << " mm" << std::endl << std::endl
              << "Direction:" << std::endl << std::endl
              << "v.x = "  << v.x() << std::endl
              << "v.y = "  << v.y() << std::endl
              << "v.z = "  << v.z() << std::endl << std::endl
              << "Proposed distance :" << std::endl << std::endl
              << "snxt = "    << snxt << " mm" << std::endl;
      message.precision(oldprc);
      UUtils::Exception("USphere::DistanceToOut(p,v,..)",
                        "GeomSolids1002", UWarning, 1, message.str().c_str());
      break;
  }
  if (snxt == UUtils::Infinity())
  {
    cout << std::endl;

    std::ostringstream message;
    int oldprc = message.precision(16);
    message << "Logic error: snxt = UUtils::Infinity() ???" << std::endl
            << "Position:"  << std::endl << std::endl
            << "p.x = "  << p.x() << " mm" << std::endl
            << "p.y = "  << p.y() << " mm" << std::endl
            << "p.z = "  << p.z() << " mm" << std::endl << std::endl
            << "Rp = " << std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z())
            << " mm" << std::endl << std::endl
            << "Direction:" << std::endl << std::endl
            << "v.x = "  << v.x() << std::endl
            << "v.y = "  << v.y() << std::endl
            << "v.z = "  << v.z() << std::endl << std::endl
            << "Proposed distance :" << std::endl << std::endl
            << "snxt = "    << snxt << " mm" << std::endl;
    message.precision(oldprc);
    UUtils::Exception("USphere::DistanceToOut(p,v,..)",
                      "GeomSolids1002", UWarning, 1, message.str().c_str());
  }

  return snxt;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

double USphere::SafetyFromInside(const UVector3& p, bool /*aAccurate*/) const
{
  double safe = 0.0, safeRMin, safeRMax, safePhi, safeTheta;
  double rho2, rds, rho;
  double pTheta, dTheta1 = UUtils::Infinity(),dTheta2 = UUtils::Infinity();
  rho2 = p.x() * p.x() + p.y() * p.y();
  rds = std::sqrt(rho2 + p.z() * p.z());
  rho = std::sqrt(rho2);

#ifdef UDEBUG
  if (Inside(p) == eOutside)
  {
    int old_prc = cout.precision(16);
    cout << std::endl;

    cout << "Position:"  << std::endl << std::endl;
    cout << "p.x = "  << p.x() << " mm" << std::endl;
    cout << "p.y = "  << p.y() << " mm" << std::endl;
    cout << "p.z = "  << p.z() << " mm" << std::endl << std::endl;
    cout.precision(old_prc);
    UUtils::Exception("USphere::DistanceToOut(p)",
                      "GeomSolids1002", UWarning, 1, "Point p is outside !?");
  }
#endif

  // Distance to r shells
  //
  safeRMax = fRmax-rds;
  safe = safeRMax;  
  if (fRmin)
  {
     safeRMin = rds-fRmin;
     safe = std::min( safeRMin, safeRMax ); 
  }

  // Distance to phi extent
  //
  if ( !fFullPhiSphere )
  {
     if (rho>0.0)
     {
        if ((p.y()*cosCPhi-p.x()*sinCPhi)<=0)
        {
           safePhi=-(p.x()*sinSPhi-p.y()*cosSPhi);
        }
        else
        {
           safePhi=(p.x()*sinEPhi-p.y()*cosEPhi);
        }
     }
     else
     {
        safePhi = 0.0;  // Distance to both Phi surfaces (extended)
     }
     // Both cases above can be improved - in case fRMin > 0.0
     //  although it may be costlier (good for precise, not fast version)
     
     safe= std::min(safe, safePhi);
  }

  // Distance to Theta extent
  //
  if ( !fFullThetaSphere )
  {
    if( rds > 0.0 )
    {
       pTheta=std::acos(p.z()/rds);
       if (pTheta<0) { pTheta+=UUtils::kPi; }
       if(fSTheta>0.)
       { dTheta1=pTheta-fSTheta;}
       if(eTheta<UUtils::kPi)
       { dTheta2=eTheta-pTheta;}
      
       safeTheta=rds*std::sin(std::min(dTheta1, dTheta2) );
    }
    else
    {
       safeTheta= 0.0;
         // An improvement will be to return negative answer if outside (TODO)
    }
    safe = std::min( safe, safeTheta );
  }

  if (safe<0.0) { safe=0; }
    // An improvement to return negative answer if outside (TODO)
  
  return safe;
}

//////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz Cross section
//          [4-7] +fDz Cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility
//  Potential improvement: For last slice, use actual ending angle
//                         to avoid rounding error problems.

/*
UVector3List*
USphere::CreateRotatedVertices(const UAffineTransform& pTransform,
                                       int& noPolygonVertices) const
{
  UVector3List *vertices;
  UVector3 vertex;
  double meshAnglePhi,meshRMax,crossAnglePhi,
           coscrossAnglePhi,sincrossAnglePhi,sAnglePhi;
  double meshTheta,crossTheta,startTheta;
  double rMaxX,rMaxY,rMinX,rMinY,rMinZ,rMaxZ;
  int crossSectionPhi,noPhiCrossSections,crossSectionTheta,noThetaSections;

  // Phi Cross sections

  noPhiCrossSections = int(fDPhi/UUtils::kMeshAngleDefault)+1;

  if (noPhiCrossSections<UUtils::kMinMeshSections)
  {
    noPhiCrossSections=UUtils::kMinMeshSections;
  }
  else if (noPhiCrossSections>UUtils::kMaxMeshSections)
  {
    noPhiCrossSections=UUtils::kMaxMeshSections;
  }
  meshAnglePhi=fDPhi/(noPhiCrossSections-1);

  // If complete in phi, Set start angle such that mesh will be at fRMax
  // on the x axis. Will give better extent calculations when not rotated.

  if (fFullPhiSphere)
  {
    sAnglePhi = -meshAnglePhi*0.5;
  }
    else
  {
    sAnglePhi=fSPhi;
  }

  // Theta Cross sections

  noThetaSections = int(fDTheta/UUtils::kMeshAngleDefault)+1;

  if (noThetaSections<UUtils::kMinMeshSections)
  {
    noThetaSections=UUtils::kMinMeshSections;
  }
  else if (noThetaSections>UUtils::kMaxMeshSections)
  {
    noThetaSections=UUtils::kMaxMeshSections;
  }
  meshTheta=fDTheta/(noThetaSections-1);

  // If complete in Theta, Set start angle such that mesh will be at fRMax
  // on the z axis. Will give better extent calculations when not rotated.

  if (fFullThetaSphere)
  {
    startTheta = -meshTheta*0.5;
  }
  else
  {
    startTheta=fSTheta;
  }

  meshRMax = (meshAnglePhi >= meshTheta) ?
             fRmax/std::cos(meshAnglePhi*0.5) : fRmax/std::cos(meshTheta*0.5);
  double* cosCrossTheta = new double[noThetaSections];
  double* sinCrossTheta = new double[noThetaSections];
  vertices=new UVector3List();
  if (vertices && cosCrossTheta && sinCrossTheta)
  {
    vertices->reserve(noPhiCrossSections*(noThetaSections*2));
    for (crossSectionPhi=0;
         crossSectionPhi<noPhiCrossSections; crossSectionPhi++)
    {
      crossAnglePhi=sAnglePhi+crossSectionPhi*meshAnglePhi;
      coscrossAnglePhi=std::cos(crossAnglePhi);
      sincrossAnglePhi=std::sin(crossAnglePhi);
      for (crossSectionTheta=0;
           crossSectionTheta<noThetaSections;crossSectionTheta++)
      {
        // Compute coordinates of Cross section at section crossSectionPhi
        //
        crossTheta=startTheta+crossSectionTheta*meshTheta;
        cosCrossTheta[crossSectionTheta]=std::cos(crossTheta);
        sinCrossTheta[crossSectionTheta]=std::sin(crossTheta);

        rMinX=fRmin*sinCrossTheta[crossSectionTheta]*coscrossAnglePhi;
        rMinY=fRmin*sinCrossTheta[crossSectionTheta]*sincrossAnglePhi;
        rMinZ=fRmin*cosCrossTheta[crossSectionTheta];

        vertex=UVector3(rMinX,rMinY,rMinZ);
        vertices->push_back(pTransform.TransformPoint(vertex));

      }   // Theta forward

      for (crossSectionTheta=noThetaSections-1;
           crossSectionTheta>=0; crossSectionTheta--)
      {
        rMaxX=meshRMax*sinCrossTheta[crossSectionTheta]*coscrossAnglePhi;
        rMaxY=meshRMax*sinCrossTheta[crossSectionTheta]*sincrossAnglePhi;
        rMaxZ=meshRMax*cosCrossTheta[crossSectionTheta];

        vertex=UVector3(rMaxX,rMaxY,rMaxZ);
        vertices->push_back(pTransform.TransformPoint(vertex));

      }  // Theta back
    }  // Phi
    noPolygonVertices = noThetaSections*2;
  }
  else
  {

  UUtils::Exception("USphere::CreateRotatedVertices()",
                "GeomSolids0003", UFatalError,1,
                "Error in allocation of vertices. Out of memory !");
  }

  delete [] cosCrossTheta;
  delete [] sinCrossTheta;

  return vertices;
}
*/

//////////////////////////////////////////////////////////////////////////
//
// UEntityType

UGeometryType USphere::GetEntityType() const
{
  return std::string("Sphere");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
VUSolid* USphere::Clone() const
{
  return new USphere(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& USphere::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: USphere\n"
     << " Parameters: \n"
     << "                inner radius: " << fRmin << " mm \n"
     << "                outer radius: " << fRmax << " mm \n"
     << "                starting phi of segment  : " << fSPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "                delta phi of segment     : " << fDPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "                starting theta of segment: " << fSTheta / (UUtils::kPi / 180.0) << " degrees \n"
     << "                delta theta of segment   : " << fDTheta / (UUtils::kPi / 180.0) << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

UVector3 USphere::GetPointOnSurface() const
{
  double zRand, aOne, aTwo, aThr, aFou, aFiv, chose, phi, sinphi, cosphi;
  double height1, height2, slant1, slant2, costheta, sintheta, rRand;

  height1 = (fRmax - fRmin) * cosSTheta;
  height2 = (fRmax - fRmin) * cosETheta;
  slant1  = std::sqrt(UUtils::sqr((fRmax - fRmin) * sinSTheta) + height1 * height1);
  slant2  = std::sqrt(UUtils::sqr((fRmax - fRmin) * sinETheta) + height2 * height2);
  rRand  = UUtils::GetRadiusInRing(fRmin, fRmax);

  aOne = fRmax * fRmax * fDPhi * (cosSTheta - cosETheta);
  aTwo = fRmin * fRmin * fDPhi * (cosSTheta - cosETheta);
  aThr = fDPhi * ((fRmax + fRmin) * sinSTheta) * slant1;
  aFou = fDPhi * ((fRmax + fRmin) * sinETheta) * slant2;
  aFiv = 0.5 * fDTheta * (fRmax * fRmax - fRmin * fRmin);

  phi = UUtils::Random(fSPhi, ePhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);
  costheta = UUtils::Random(cosETheta, cosSTheta);
  sintheta = std::sqrt(1. - UUtils::sqr(costheta));

  if (fFullPhiSphere)
  {
    aFiv = 0;
  }
  if (fSTheta == 0)
  {
    aThr = 0;
  }
  if (eTheta == UUtils::kPi)
  {
    aFou = 0;
  }
  if (fSTheta == UUtils::kPi / 2)
  {
    aThr = UUtils::kPi * (fRmax * fRmax - fRmin * fRmin);
  }
  if (eTheta == UUtils::kPi / 2)
  {
    aFou = UUtils::kPi * (fRmax * fRmax - fRmin * fRmin);
  }

  chose = UUtils::Random(0., aOne + aTwo + aThr + aFou + 2.*aFiv);
  if ((chose >= 0.) && (chose < aOne))
  {
    return UVector3(fRmax * sintheta * cosphi,
                    fRmax * sintheta * sinphi, fRmax * costheta);
  }
  else if ((chose >= aOne) && (chose < aOne + aTwo))
  {
    return UVector3(fRmin * sintheta * cosphi,
                    fRmin * sintheta * sinphi, fRmin * costheta);
  }
  else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aThr))
  {
    if (fSTheta != UUtils::kPi / 2)
    {
      zRand = UUtils::Random(fRmin * cosSTheta, fRmax * cosSTheta);
      return UVector3(tanSTheta * zRand * cosphi,
                      tanSTheta * zRand * sinphi, zRand);
    }
    else
    {
      return UVector3(rRand * cosphi, rRand * sinphi, 0.);
    }
  }
  else if ((chose >= aOne + aTwo + aThr) && (chose < aOne + aTwo + aThr + aFou))
  {
    if (eTheta != UUtils::kPi / 2)
    {
      zRand = UUtils::Random(fRmin * cosETheta, fRmax * cosETheta);
      return UVector3(tanETheta * zRand * cosphi,
                      tanETheta * zRand * sinphi, zRand);
    }
    else
    {
      return UVector3(rRand * cosphi, rRand * sinphi, 0.);
    }
  }
  else if ((chose >= aOne + aTwo + aThr + aFou) && (chose < aOne + aTwo + aThr + aFou + aFiv))
  {
    return UVector3(rRand * sintheta * cosSPhi,
                    rRand * sintheta * sinSPhi, rRand * costheta);
  }
  else
  {
    return UVector3(rRand * sintheta * cosEPhi,
                    rRand * sintheta * sinEPhi, rRand * costheta);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// GetSurfaceArea

double USphere::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    double Rsq = fRmax * fRmax;
    double rsq = fRmin * fRmin;

    fSurfaceArea = fDPhi * (rsq + Rsq) * (cosSTheta - cosETheta);
    if (!fFullPhiSphere)
    {
      fSurfaceArea = fSurfaceArea + fDTheta * (Rsq - rsq);
    }
    if (fSTheta > 0)
    {
      double acos1 = std::acos(std::pow(sinSTheta, 2) * std::cos(fDPhi)
                               + std::pow(cosSTheta, 2));
      if (fDPhi > UUtils::kPi)
      {
        fSurfaceArea = fSurfaceArea + 0.5 * (Rsq - rsq) * (2 * UUtils::kPi - acos1);
      }
      else
      {
        fSurfaceArea = fSurfaceArea + 0.5 * (Rsq - rsq) * acos1;
      }
    }
    if (eTheta < UUtils::kPi)
    {
      double acos2 = std::acos(std::pow(sinETheta, 2) * std::cos(fDPhi)
                               + std::pow(cosETheta, 2));
      if (fDPhi > UUtils::kPi)
      {
        fSurfaceArea = fSurfaceArea + 0.5 * (Rsq - rsq) * (2 * UUtils::kPi - acos2);
      }
      else
      {
        fSurfaceArea = fSurfaceArea + 0.5 * (Rsq - rsq) * acos2;
      }
    }
  }
  return fSurfaceArea;
}



void USphere::Extent(UVector3& aMin, UVector3& aMax) const
{
  aMin.Set(-fRmax);
  aMax.Set(fRmax);
}

void USphere::GetParametersList(int, double* aArray) const
{
  aArray[0] = GetInnerRadius();
  aArray[1] = GetOuterRadius();
  aArray[2] = GetStartPhiAngle();
  aArray[3] = GetDeltaPhiAngle();
  aArray[4] = GetStartThetaAngle();
  aArray[5] = GetDeltaThetaAngle();
}
