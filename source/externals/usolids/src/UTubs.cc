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
// UTubs
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UTubs.hh"

using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

UTubs::UTubs(const std::string& pName,
             double pRMin, double pRMax,
             double pDz,
             double pSPhi, double pDPhi)
  : VUSolid(pName.c_str()), fRMin(pRMin), fRMax(pRMax), fDz(pDz), fSPhi(0), fDPhi(0)
{
  fCubicVolume=0.;
  fSurfaceArea=0.;
  kRadTolerance = frTolerance;
  kAngTolerance = faTolerance;

  if (pDz <= 0) // Check z-len
  {
    std::ostringstream message;
    message << "Negative Z half-length (" << pDz << ") in solid: " << GetName();
    UUtils::Exception("UTubs::UTubs()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  if ((pRMin >= pRMax) || (pRMin < 0))   // Check radii
  {
    std::ostringstream message;
    message << "Invalid values for radii in solid: " << GetName()
            << std::endl
            << "pRMin = " << pRMin << ", pRMax = " << pRMax;
    UUtils::Exception("UTubs::UTubs()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Check angles
  //
  CheckPhiAngles(pSPhi, pDPhi);
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
UTubs::UTubs()
  : VUSolid(""), fCubicVolume(0.), fSurfaceArea(0.),
    kRadTolerance(0.), kAngTolerance(0.),
    fRMin(0.), fRMax(0.), fDz(0.), fSPhi(0.), fDPhi(0.),
    fSinCPhi(0.), fCosCPhi(0.), fCosHDPhiOT(0.), fCosHDPhiIT(0.),
    fSinSPhi(0.), fCosSPhi(0.), fSinEPhi(0.), fCosEPhi(0.),
    fSinSPhiDPhi(0.), fCosSPhiDPhi(0.),
    fPhiFullTube(false)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

UTubs::~UTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

UTubs::UTubs(const UTubs& rhs)
  : VUSolid(rhs),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    kRadTolerance(rhs.kRadTolerance), kAngTolerance(rhs.kAngTolerance),
    fRMin(rhs.fRMin), fRMax(rhs.fRMax), fDz(rhs.fDz),
    fSPhi(rhs.fSPhi), fDPhi(rhs.fDPhi),
    fSinCPhi(rhs.fSinCPhi), fCosCPhi(rhs.fCosCPhi),
    fCosHDPhiOT(rhs.fCosHDPhiOT), fCosHDPhiIT(rhs.fCosHDPhiIT),
    fSinSPhi(rhs.fSinSPhi), fCosSPhi(rhs.fCosSPhi),
    fSinEPhi(rhs.fSinEPhi), fCosEPhi(rhs.fCosEPhi),
    fSinSPhiDPhi(rhs.fSinSPhiDPhi), fCosSPhiDPhi(rhs.fCosSPhiDPhi),
    fPhiFullTube(rhs.fPhiFullTube)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

UTubs& UTubs::operator = (const UTubs& rhs)
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
  kRadTolerance = rhs.kRadTolerance;
  kAngTolerance = rhs.kAngTolerance;
  fRMin = rhs.fRMin;
  fRMax = rhs.fRMax;
  fDz = rhs.fDz;
  fSPhi = rhs.fSPhi;
  fDPhi = rhs.fDPhi;
  fSinCPhi = rhs.fSinCPhi;
  fCosCPhi = rhs.fCosCPhi;
  fCosHDPhiOT = rhs.fCosHDPhiOT;
  fCosHDPhiIT = rhs.fCosHDPhiIT;
  fSinSPhi = rhs.fSinSPhi;
  fCosSPhi = rhs.fCosSPhi;
  fSinEPhi = rhs.fSinEPhi;
  fCosEPhi = rhs.fCosEPhi;
  fSinSPhiDPhi = rhs.fSinSPhiDPhi;
  fCosSPhiDPhi = rhs.fCosSPhiDPhi;

  fPhiFullTube = rhs.fPhiFullTube;

  return *this;
}


///////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

inline VUSolid::EnumInside UTubs::Inside(const UVector3& p) const
{
  double r2, pPhi, tolRMin, tolRMax;
  VUSolid::EnumInside in = eOutside;
  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfRadTolerance = kRadTolerance * 0.5;
  static const double halfAngTolerance = kAngTolerance * 0.5;

  if (std::fabs(p.z()) <= fDz - halfCarTolerance)
  {
    r2 = p.x() * p.x() + p.y() * p.y();

    if (fRMin)
    {
      tolRMin = fRMin + halfRadTolerance;
    }
    else
    {
      tolRMin = 0;
    }

    tolRMax = fRMax - halfRadTolerance;

    if ((r2 >= tolRMin * tolRMin) && (r2 <= tolRMax * tolRMax))
    {
      if (fPhiFullTube)
      {
        in = eInside;
      }
      else
      {
        // Try inner tolerant phi boundaries (=>inside)
        // if not inside, try outer tolerant phi boundaries

        if ((tolRMin == 0) && (std::fabs(p.x()) <= halfCarTolerance)
            && (std::fabs(p.y()) <= halfCarTolerance))
        {
          in = eSurface;
        }
        else
        {
          pPhi = std::atan2(p.y(), p.x());
          if (pPhi < -halfAngTolerance)
          {
            pPhi += 2 * UUtils::kPi;  // 0<=pPhi<2UUtils::kPi
          }

          if (fSPhi >= 0)
          {
            if ((std::fabs(pPhi) < halfAngTolerance)
                && (std::fabs(fSPhi + fDPhi - 2 * UUtils::kPi) < halfAngTolerance))
            {
              pPhi += 2 * UUtils::kPi; // 0 <= pPhi < 2UUtils::kPi
            }
            if ((pPhi >= fSPhi + halfAngTolerance)
                && (pPhi <= fSPhi + fDPhi - halfAngTolerance))
            {
              in = eInside;
            }
            else if ((pPhi >= fSPhi - halfAngTolerance)
                     && (pPhi <= fSPhi + fDPhi + halfAngTolerance))
            {
              in = eSurface;
            }
          }
          else  // fSPhi < 0
          {
            if ((pPhi <= fSPhi + 2 * UUtils::kPi - halfAngTolerance)
                && (pPhi >= fSPhi + fDPhi + halfAngTolerance))
            {
              ; //eOutside
            }
            else if ((pPhi <= fSPhi + 2 * UUtils::kPi + halfAngTolerance)
                     && (pPhi >= fSPhi + fDPhi - halfAngTolerance))
            {
              in = eSurface;
            }
            else
            {
              in = eInside;
            }
          }
        }
      }
    }
    else  // Try generous boundaries
    {
      tolRMin = fRMin - halfRadTolerance;
      tolRMax = fRMax + halfRadTolerance;

      if (tolRMin < 0)
      {
        tolRMin = 0;
      }

      if ((r2 >= tolRMin * tolRMin) && (r2 <= tolRMax * tolRMax))
      {
        if (fPhiFullTube || (r2 <= halfRadTolerance * halfRadTolerance))
        {
          // Continuous in phi or on z-axis
          in = eSurface;
        }
        else // Try outer tolerant phi boundaries only
        {
          pPhi = std::atan2(p.y(), p.x());

          if (pPhi < -halfAngTolerance)
          {
            pPhi += 2 * UUtils::kPi;  // 0<=pPhi<2UUtils::kPi
          }
          if (fSPhi >= 0)
          {
            if ((std::fabs(pPhi) < halfAngTolerance)
                && (std::fabs(fSPhi + fDPhi - 2 * UUtils::kPi) < halfAngTolerance))
            {
              pPhi += 2 * UUtils::kPi; // 0 <= pPhi < 2UUtils::kPi
            }
            if ((pPhi >= fSPhi - halfAngTolerance)
                && (pPhi <= fSPhi + fDPhi + halfAngTolerance))
            {
              in = eSurface;
            }
          }
          else  // fSPhi < 0
          {
            if ((pPhi <= fSPhi + 2 * UUtils::kPi - halfAngTolerance)
                && (pPhi >= fSPhi + fDPhi + halfAngTolerance))
            {
              ; // eOutside
            }
            else
            {
              in = eSurface;
            }
          }
        }
      }
    }
  }
  else if (std::fabs(p.z()) <= fDz + halfCarTolerance)
  {
    // Check within tolerant r limits
    r2 = p.x() * p.x() + p.y() * p.y();
    tolRMin = fRMin - halfRadTolerance;
    tolRMax = fRMax + halfRadTolerance;

    if (tolRMin < 0)
    {
      tolRMin = 0;
    }

    if ((r2 >= tolRMin * tolRMin) && (r2 <= tolRMax * tolRMax))
    {
      if (fPhiFullTube || (r2 <= halfRadTolerance * halfRadTolerance))
      {
        // Continuous in phi or on z-axis
        in = eSurface;
      }
      else // Try outer tolerant phi boundaries
      {
        pPhi = std::atan2(p.y(), p.x());

        if (pPhi < -halfAngTolerance)
        {
          pPhi += 2 * UUtils::kPi;  // 0<=pPhi<2UUtils::kPi
        }
        if (fSPhi >= 0)
        {
          if ((std::fabs(pPhi) < halfAngTolerance)
              && (std::fabs(fSPhi + fDPhi - 2 * UUtils::kPi) < halfAngTolerance))
          {
            pPhi += 2 * UUtils::kPi; // 0 <= pPhi < 2UUtils::kPi
          }
          if ((pPhi >= fSPhi - halfAngTolerance)
              && (pPhi <= fSPhi + fDPhi + halfAngTolerance))
          {
            in = eSurface;
          }
        }
        else  // fSPhi < 0
        {
          if ((pPhi <= fSPhi + 2 * UUtils::kPi - halfAngTolerance)
              && (pPhi >= fSPhi + fDPhi + halfAngTolerance))
          {
            ;
          }
          else
          {
            in = eSurface;
          }
        }
      }
    }
  }
  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return Unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

bool UTubs::Normal(const UVector3& p, UVector3& n) const
{
  int noSurfaces = 0;
  double rho, pPhi;
  double distZ, distRMin, distRMax;
  double distSPhi = UUtils::kInfinity, distEPhi = UUtils::kInfinity;

  static const double halfCarTolerance = 0.5 * VUSolid::Tolerance();
  static const double halfAngTolerance = 0.5 * kAngTolerance;

  UVector3 norm, sumnorm(0., 0., 0.);
  UVector3 nZ = UVector3(0, 0, 1.0);
  UVector3 nR, nPs, nPe;

  rho = std::sqrt(p.x() * p.x() + p.y() * p.y());

  distRMin = std::fabs(rho - fRMin);
  distRMax = std::fabs(rho - fRMax);
  distZ   = std::fabs(std::fabs(p.z()) - fDz);

  if (!fPhiFullTube)    // Protected against (0,0,z)
  {
    if (rho > halfCarTolerance)
    {
      pPhi = std::atan2(p.y(), p.x());

      if (pPhi < fSPhi - halfCarTolerance)
      {
        pPhi += 2 * UUtils::kPi;
      }
      else if (pPhi > fSPhi + fDPhi + halfCarTolerance)
      {
        pPhi -= 2 * UUtils::kPi;
      }

      distSPhi = std::fabs(pPhi - fSPhi);
      distEPhi = std::fabs(pPhi - fSPhi - fDPhi);
    }
    else if (!fRMin)
    {
      distSPhi = 0.;
      distEPhi = 0.;
    }
    nPs = UVector3(fSinSPhi, -fCosSPhi, 0);
    nPe = UVector3(-fSinSPhiDPhi /* std::sin(fSPhi+fDPhi)*/, fCosSPhiDPhi, 0);
  }
  if (rho > halfCarTolerance)
  {
    nR = UVector3(p.x() / rho, p.y() / rho, 0);
  }

  if (distRMax <= halfCarTolerance)
  {
    noSurfaces ++;
    sumnorm += nR;
  }
  if (fRMin && (distRMin <= halfCarTolerance))
  {
    noSurfaces ++;
    sumnorm -= nR;
  }
  if (fDPhi < 2 * UUtils::kPi)
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
  if (distZ <= halfCarTolerance)
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
    UUtils::Exception("UTubs::SurfaceNormal(p)", "GeomSolids1002",
                      UWarning, 1, "Point p is not on surface !?");
    int oldprc = cout.precision(20);
    cout << "UTubs::SN ( " << p.x() << ", " << p.y() << ", " << p.z() << " ); "
         << std::endl << std::endl;
    cout.precision(oldprc);
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

  return noSurfaces; // TODO: return true or false on validity
}

/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

UVector3 UTubs::ApproxSurfaceNormal(const UVector3& p) const
{
  ENorm side;
  UVector3 norm;
  double rho, phi;
  double distZ, distRMin, distRMax, distSPhi, distEPhi, distMin;

  rho = std::sqrt(p.x() * p.x() + p.y() * p.y());

  distRMin = std::fabs(rho - fRMin);
  distRMax = std::fabs(rho - fRMax);
  distZ   = std::fabs(std::fabs(p.z()) - fDz);

  if (distRMin < distRMax) // First minimum
  {
    if (distZ < distRMin)
    {
      distMin = distZ;
      side    = kNZ;
    }
    else
    {
      distMin = distRMin;
      side    = kNRMin  ;
    }
  }
  else
  {
    if (distZ < distRMax)
    {
      distMin = distZ;
      side    = kNZ ;
    }
    else
    {
      distMin = distRMax;
      side    = kNRMax  ;
    }
  }
  if (!fPhiFullTube &&  rho)  // Protected against (0,0,z)
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

    if (distSPhi < distEPhi) // Find new minimum
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
    case kNRMin : // Inner radius
    {
      norm = UVector3(-p.x() / rho, -p.y() / rho, 0);
      break;
    }
    case kNRMax : // Outer radius
    {
      norm = UVector3(p.x() / rho, p.y() / rho, 0);
      break;
    }
    case kNZ :    // + or - dz
    {
      if (p.z() > 0)
      {
        norm = UVector3(0, 0, 1);
      }
      else
      {
        norm = UVector3(0, 0, -1);
      }
      break;
    }
    case kNSPhi:
    {
      norm = UVector3(fSinSPhi, -fCosSPhi, 0);
      break;
    }
    case kNEPhi:
    {
      norm = UVector3(-fSinSPhiDPhi, fCosSPhiDPhi, 0);
      break;
    }
    default:      // Should never reach this case ...
    {
      // DumpInfo();
      UUtils::Exception("UTubs::ApproxSurfaceNormal()",
                        "GeomSolids1002", UWarning, 1,
                        "Undefined side for valid surface normal to solid.");
      break;
    }
  }
  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return UUtils::kInfinity if no intersection, or intersection distance <= tolerance
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
// - 'if valid' implies tolerant checking of intersection points

double UTubs::DistanceToIn(const UVector3& p, const UVector3& v, double) const
{
  double snxt = UUtils::kInfinity;      // snxt = default return value
  double tolORMin2, tolIRMax2;  // 'generous' radii squared
  double tolORMax2, tolIRMin2, tolODz, tolIDz;
  const double dRmax = 100.*fRMax;

  static const double halfCarTolerance = 0.5 * VUSolid::Tolerance();
  static const double halfRadTolerance = 0.5 * kRadTolerance;

  // Intersection point variables
  //
  double Dist, sd, xi, yi, zi, rho2, inum, iden, cosPsi, Comp;
  double t1, t2, t3, b, c, d;    // Quadratic solver variables

  // Calculate tolerant rmin and rmax

  if (fRMin > kRadTolerance)
  {
    tolORMin2 = (fRMin - halfRadTolerance) * (fRMin - halfRadTolerance);
    tolIRMin2 = (fRMin + halfRadTolerance) * (fRMin + halfRadTolerance);
  }
  else
  {
    tolORMin2 = 0.0;
    tolIRMin2 = 0.0;
  }
  tolORMax2 = (fRMax + halfRadTolerance) * (fRMax + halfRadTolerance);
  tolIRMax2 = (fRMax - halfRadTolerance) * (fRMax - halfRadTolerance);

  // Intersection with Z surfaces

  tolIDz = fDz - halfCarTolerance;
  tolODz = fDz + halfCarTolerance;

  if (std::fabs(p.z()) >= tolIDz)
  {
    if (p.z() * v.z() < 0)   // at +Z going in -Z or visa versa
    {
      sd = (std::fabs(p.z()) - fDz) / std::fabs(v.z()); // Z intersect distance

      if (sd < 0.0)
      {
        sd = 0.0;
      }

      xi   = p.x() + sd * v.x();              // Intersection coords
      yi   = p.y() + sd * v.y();
      rho2 = xi * xi + yi * yi;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi * fCosCPhi + yi * fSinCPhi;
          iden   = std::sqrt(rho2);
          cosPsi = inum / iden;
          if (cosPsi >= fCosHDPhiIT)
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
    else
    {
      if (snxt < halfCarTolerance)
      {
        snxt = 0;
      }
      return snxt;  // On/outside extent, and heading away
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

  t1 = 1.0 - v.z() * v.z();
  t2 = p.x() * v.x() + p.y() * v.y();
  t3 = p.x() * p.x() + p.y() * p.y();

  if (t1 > 0)          // Check not || to z axis
  {
    b = t2 / t1;
    c = t3 - fRMax * fRMax;
    if ((t3 >= tolORMax2) && (t2 < 0)) // This also handles the tangent case
    {
      // Try outer cylinder intersection
      //          c=(t3-fRMax*fRMax)/t1;

      c /= t1;
      d = b * b - c;

      if (d >= 0) // If real root
      {
        sd = c / (-b + std::sqrt(d));
        if (sd >= 0)  // If 'forwards'
        {
          if (sd > dRmax) // Avoid rounding errors due to precision issues on
          {
            // 64 bits systems. Split long distances and recompute
            double fTerm = sd - std::fmod(sd, dRmax);
            sd = fTerm + DistanceToIn(p + fTerm * v, v);
          }
          // Check z intersection
          //
          zi = p.z() + sd * v.z();
          if (std::fabs(zi) <= tolODz)
          {
            // Z ok. Check phi intersection if reqd
            //
            if (fPhiFullTube)
            {
              return sd;
            }
            else
            {
              xi     = p.x() + sd * v.x();
              yi     = p.y() + sd * v.y();
              cosPsi = (xi * fCosCPhi + yi * fSinCPhi) / fRMax;
              if (cosPsi >= fCosHDPhiIT)
              {
                return sd;
              }
            }
          } //  end if std::fabs(zi)
        }   //  end if (sd>=0)
      }     //  end if (d>=0)
    }       //  end if (r>=fRMax)
    else
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)

      if ((t3 > tolIRMin2) && (t2 < 0) && (std::fabs(p.z()) <= tolIDz))
      {
        // Inside both radii, delta r -ve, inside z extent

        if (!fPhiFullTube)
        {
          inum   = p.x() * fCosCPhi + p.y() * fSinCPhi;
          iden   = std::sqrt(t3);
          cosPsi = inum / iden;
          if (cosPsi >= fCosHDPhiIT)
          {
            // In the old version, the small negative tangent for the point
            // on surface was not taken in account, and returning 0.0 ...
            // New version: check the tangent for the point on surface and
            // if no intersection, return UUtils::kInfinity, if intersection instead
            // return sd.
            //
            c = t3 - fRMax * fRMax;
            if (c <= 0.0)
            {
              return 0.0;
            }
            else
            {
              c = c / t1;
              d = b * b - c;
              if (d >= 0.0)
              {
                snxt = c / (-b + std::sqrt(d)); // using safe solution
                // for quadratic equation
                if (snxt < halfCarTolerance)
                {
                  snxt = 0;
                }
                return snxt;
              }
              else
              {
                return UUtils::kInfinity;
              }
            }
          }
        }
        else
        {
          // In the old version, the small negative tangent for the point
          // on surface was not taken in account, and returning 0.0 ...
          // New version: check the tangent for the point on surface and
          // if no intersection, return UUtils::kInfinity, if intersection instead
          // return sd.
          //
          c = t3 - fRMax * fRMax;
          if (c <= 0.0)
          {
            return 0.0;
          }
          else
          {
            c = c / t1;
            d = b * b - c;
            if (d >= 0.0)
            {
              snxt = c / (-b + std::sqrt(d)); // using safe solution
              // for quadratic equation
              if (snxt < halfCarTolerance)
              {
                snxt = 0;
              }
              return snxt;
            }
            else
            {
              return UUtils::kInfinity;
            }
          }
        } // end if  (!fPhiFullTube)
      }  // end if   (t3>tolIRMin2)
    }    // end if   (Inside Outer Radius)
    if (fRMin)     // Try inner cylinder intersection
    {
      c = (t3 - fRMin * fRMin) / t1;
      d = b * b - c;
      if (d >= 0.0)    // If real root
      {
        // Always want 2nd root - we are outside and know rmax Hit was bad
        // - If on surface of rmin also need farthest root

        sd = (b > 0.) ? c / (-b - std::sqrt(d)) : (-b + std::sqrt(d));
        if (sd >= -halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if (sd < 0.0)
          {
            sd = 0.0;
          }
          if (sd > dRmax) // Avoid rounding errors due to precision issues seen
          {
            // 64 bits systems. Split long distances and recompute
            double fTerm = sd - std::fmod(sd, dRmax);
            sd = fTerm + DistanceToIn(p + fTerm * v, v);
          }
          zi = p.z() + sd * v.z();
          if (std::fabs(zi) <= tolODz)
          {
            // Z ok. Check phi
            //
            if (fPhiFullTube)
            {
              return sd;
            }
            else
            {
              xi     = p.x() + sd * v.x();
              yi     = p.y() + sd * v.y();
              cosPsi = (xi * fCosCPhi + yi * fSinCPhi) / fRMin;
              if (cosPsi >= fCosHDPhiIT)
              {
                // Good inner radius isect
                // - but earlier phi isect still possible

                snxt = sd;
              }
            }
          }       //    end if std::fabs(zi)
        }         //    end if (sd>=0)
      }           //    end if (d>=0)
    }             //    end if (fRMin)
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to VUSolid::Tolerance()*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?
  //
  if (!fPhiFullTube)
  {
    // First phi surface (Starting phi)
    //
    Comp    = v.x() * fSinSPhi - v.y() * fCosSPhi;

    if (Comp < 0)    // Component in outwards normal dirn
    {
      Dist = (p.y() * fCosSPhi - p.x() * fSinSPhi);

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
            xi   = p.x() + sd * v.x();
            yi   = p.y() + sd * v.y();
            rho2 = xi * xi + yi * yi;

            if (((rho2 >= tolIRMin2) && (rho2 <= tolIRMax2))
                || ((rho2 > tolORMin2) && (rho2 < tolIRMin2)
                    && (v.y() * fCosSPhi - v.x() * fSinSPhi > 0)
                    && (v.x() * fCosSPhi + v.y() * fSinSPhi >= 0))
                || ((rho2 > tolIRMax2) && (rho2 < tolORMax2)
                    && (v.y() * fCosSPhi - v.x() * fSinSPhi > 0)
                    && (v.x() * fCosSPhi + v.y() * fSinSPhi < 0)))
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ((yi * fCosCPhi - xi * fSinCPhi) <= halfCarTolerance)
              {
                snxt = sd;
              }
            }
          }
        }
      }
    }

    // Second phi surface (Ending phi)

    Comp    = -(v.x() * fSinEPhi - v.y() * fCosEPhi);

    if (Comp < 0)   // Component in outwards normal dirn
    {
      Dist = -(p.y() * fCosEPhi - p.x() * fSinEPhi);

      if (Dist < halfCarTolerance)
      {
        sd = Dist / Comp;

        if (sd < snxt)
        {
          if (sd < 0)
          {
            sd = 0;
          }
          zi = p.z() + sd * v.z();
          if (std::fabs(zi) <= tolODz)
          {
            xi   = p.x() + sd * v.x();
            yi   = p.y() + sd * v.y();
            rho2 = xi * xi + yi * yi;
            if (((rho2 >= tolIRMin2) && (rho2 <= tolIRMax2))
                || ((rho2 > tolORMin2)  && (rho2 < tolIRMin2)
                    && (v.x() * fSinEPhi - v.y() * fCosEPhi > 0)
                    && (v.x() * fCosEPhi + v.y() * fSinEPhi >= 0))
                || ((rho2 > tolIRMax2) && (rho2 < tolORMax2)
                    && (v.x() * fSinEPhi - v.y() * fCosEPhi > 0)
                    && (v.x() * fCosEPhi + v.y() * fSinEPhi < 0)))
            {
              // z and r intersections good
              // - check intersecting with correct half-plane
              //
              if ((yi * fCosCPhi - xi * fSinCPhi) >= 0)
              {
                snxt = sd;
              }
            }                        //?? >=-halfCarTolerance
          }
        }
      }
    }        // Comp < 0
  }          // !fPhiFullTube
  if (snxt < halfCarTolerance)
  {
    snxt = 0;
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return UUtils::kInfinity if no intersection, or intersection distance <= tolerance
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
//   Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

double UTubs::SafetyFromOutside(const UVector3& p, bool aAccurate) const
{
  double safe = 0.0, rho, safe1, safe2, safe3;
  double safePhi;
  bool outside;

  rho  = std::sqrt(p.x() * p.x() + p.y() * p.y());
  safe1 = fRMin - rho;
  safe2 = rho - fRMax;
  safe3 = std::fabs(p.z()) - fDz;

  if (safe1 > safe2)
  {
    safe = safe1;
  }
  else
  {
    safe = safe2;
  }
  if (safe3 > safe)
  {
    safe = safe3;
  }

  if ((!fPhiFullTube) && (rho))
  {
    safePhi = SafetyToPhi(p,rho,outside);
    if ((outside) && (safePhi > safe))
    {
      safe = safePhi;
    }
  }

  if (safe < 0)
  {
    safe = 0; return safe; // point is Inside;
  }
  if (!aAccurate) return safe;
  double safsq = 0.0;
  int count = 0;
  if (safe1 > 0)
  {
    safsq += safe1 * safe1;
    count++;
  }
  if (safe2 > 0)
  {
    safsq += safe2 * safe2;
    count++;
  }
  if (safe3 > 0)
  {
    safsq += safe3 * safe3;
    count++;
  }
  if (count == 1) return safe;
  return std::sqrt(safsq);
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

// double UTubs::DistanceToOut( const UVector3& p, const UVector3& v, const bool calcNorm, bool *validNorm, UVector3 *n   ) const
double UTubs::DistanceToOut(const UVector3& p, const UVector3& v, UVector3& n, bool& validNorm, double) const
{
  ESide side = kNull , sider = kNull, sidephi = kNull;
  double snxt, srd = UUtils::kInfinity, sphi = UUtils::kInfinity, pdist;
  double deltaR, t1, t2, t3, b, c, d2, roMin2;

  static const double halfCarTolerance = VUSolid::Tolerance() * 0.5;
  static const double halfAngTolerance = kAngTolerance * 0.5;

  // Vars for phi intersection:

  double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi, roi2;

  // Z plane intersection

  if (v.z() > 0)
  {
    pdist = fDz - p.z();
    if (pdist > halfCarTolerance)
    {
      snxt = pdist / v.z();
      side = kPZ;
    }
    else
    {
      n = UVector3(0, 0, 1);
      validNorm = true;
      return snxt = 0;
    }
  }
  else if (v.z() < 0)
  {
    pdist = fDz + p.z();

    if (pdist > halfCarTolerance)
    {
      snxt = -pdist / v.z();
      side = kMZ;
    }
    else
    {
      n = UVector3(0, 0, -1);
      validNorm = true;
      return snxt = 0.0;
    }
  }
  else
  {
    snxt = UUtils::kInfinity;   // Travel perpendicular to z axis
    side = kNull;
  }

  // Radial Intersections
  //
  // Find intersection with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //
  //            t1                t2                    t3

  t1   = 1.0 - v.z() * v.z();   // since v normalised
  t2   = p.x() * v.x() + p.y() * v.y();
  t3   = p.x() * p.x() + p.y() * p.y();

  if (snxt > 10 * (fDz + fRMax))
  {
    roi2 = 2 * fRMax * fRMax;
  }
  else
  {
    roi2 = snxt * snxt * t1 + 2 * snxt * t2 + t3;  // radius^2 on +-fDz
  }

  if (t1 > 0)   // Check not parallel
  {
    // Calculate srd, r exit distance

    if ((t2 >= 0.0) && (roi2 > fRMax * (fRMax + kRadTolerance)))
    {
      // Delta r not negative => leaving via rmax

      deltaR = t3 - fRMax * fRMax;

      // NOTE: Should use rho-fRMax<-kRadTolerance*0.5
      // - avoid sqrt for efficiency

      if (deltaR < -kRadTolerance * fRMax)
      {
        b    = t2 / t1;
        c    = deltaR / t1;
        d2    = b * b - c;
        if (d2 >= 0)
        {
          srd = c / (-b - std::sqrt(d2));
        }
        else
        {
          srd = 0.;
        }
        sider = kRMax;
      }
      else
      {
        // On tolerant boundary & heading outwards (or perpendicular to)
        // outer radial surface -> leaving immediately
          n = UVector3(p.x() / fRMax, p.y() / fRMax, 0);
          validNorm = true;
          return snxt = 0; // Leaving by rmax immediately
      }
    }
    else if (t2 < 0.)   // i.e.  t2 < 0; Possible rmin intersection
    {
      roMin2 = t3 - t2 * t2 / t1; // min ro2 of the plane of movement

      if (fRMin && (roMin2 < fRMin * (fRMin - kRadTolerance)))
      {
        deltaR = t3 - fRMin * fRMin;
        b     = t2 / t1;
        c     = deltaR / t1;
        d2     = b * b - c;

        if (d2 >= 0)    // Leaving via rmin
        {
          // NOTE: SHould use rho-rmin>kRadTolerance*0.5
          // - avoid sqrt for efficiency

          if (deltaR > kRadTolerance * fRMin)
          {
            srd = c / (-b + std::sqrt(d2));
            sider = kRMin;
          }
          else
          {
            validNorm = false;  // Concave side
            n = UVector3(-p.x() / fRMin, -p.y() / fRMin, 0);
            return snxt = 0.0;
          }
        }
        else    // No rmin intersect -> must be rmax intersect
        {
          deltaR = t3 - fRMax * fRMax;
          c    = deltaR / t1;
          d2    = b * b - c;
          if (d2 >= 0.)
          {
            srd    = -b + std::sqrt(d2);
            sider = kRMax;
          }
          else // Case: On the border+t2<kRadTolerance
            //       (v is perpendicular to the surface)
          {
            n = UVector3(p.x() / fRMax, p.y() / fRMax, 0);
            validNorm = true;
            return snxt = 0.0;
          }
        }
      }
      else if (roi2 > fRMax * (fRMax + kRadTolerance))
        // No rmin intersect -> must be rmax intersect
      {
        deltaR = t3 - fRMax * fRMax;
        b     = t2 / t1;
        c     = deltaR / t1;
        d2     = b * b - c;
        if (d2 >= 0)
        {
          srd    = -b + std::sqrt(d2);
          sider = kRMax;
        }
        else // Case: On the border+t2<kRadTolerance
          //       (v is perpendicular to the surface)
        {
          n = UVector3(p.x() / fRMax, p.y() / fRMax, 0);
          validNorm = true;
          return snxt = 0.0;
        }
      }
    }

    // Phi Intersection

    if (!fPhiFullTube)
    {
      // add angle calculation with correction
      // of the difference in domain of atan2 and Sphi
      //
      vphi = std::atan2(v.y(), v.x());

      if (vphi < fSPhi - halfAngTolerance)
      {
        vphi += 2 * UUtils::kPi;
      }
      else if (vphi > fSPhi + fDPhi + halfAngTolerance)
      {
        vphi -= 2 * UUtils::kPi;
      }


      if (p.x() || p.y())    // Check if on z axis (rho not needed later)
      {
        // pDist -ve when inside

        pDistS = p.x() * fSinSPhi - p.y() * fCosSPhi;
        pDistE = -p.x() * fSinEPhi + p.y() * fCosEPhi;

        // Comp -ve when in direction of outwards normal

        compS  = -fSinSPhi * v.x() + fCosSPhi * v.y();
        compE  =  fSinEPhi * v.x() - fCosEPhi * v.y();

        sidephi = kNull;

        if (((fDPhi <= UUtils::kPi) && ((pDistS <= halfCarTolerance)
                                        && (pDistE <= halfCarTolerance)))
            || ((fDPhi >  UUtils::kPi) && !((pDistS > halfCarTolerance)
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
              if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
              {
                sidephi = kSPhi;
                if (((fSPhi - halfAngTolerance) <= vphi)
                    && ((fSPhi + fDPhi + halfAngTolerance) >= vphi))
                {
                  sphi = UUtils::kInfinity;
                }
              }
              else if (yi * fCosCPhi - xi * fSinCPhi >= 0)
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

              if ((std::fabs(xi) <= VUSolid::Tolerance()) && (std::fabs(yi) <= VUSolid::Tolerance()))
              {
                // Leaving via ending phi
                //
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
              else    // Check intersecting with correct half-plane

                if ((yi * fCosCPhi - xi * fSinCPhi) >= 0)
                {
                  // Leaving via ending phi
                  //
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
          sidephi = kSPhi; // arbitrary
          sphi    = 0.0;
        }
      }
      if (sphi < snxt)  // Order intersecttions
      {
        snxt = sphi;
        side = sidephi;
      }
    }
    if (srd < snxt) // Order intersections
    {
      snxt = srd;
      side = sider;
    }
  }
  //  if (calcNorm)
  {
    switch (side)
    {
      case kRMax:
        // Note: returned vector not normalised
        // (divide by fRMax for Unit vector)
        //
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        n = UVector3(xi / fRMax, yi / fRMax, 0);
        validNorm = true;
        break;

      case kRMin:
        xi = p.x() + snxt * v.x();
        yi = p.y() + snxt * v.y();
        n = UVector3(-xi / fRMin, -yi / fRMin, 0);
        validNorm = false;  // Rmin is inconvex
        break;

      case kSPhi:
        if (fDPhi <= UUtils::kPi)
        {
          n = UVector3(fSinSPhi, -fCosSPhi, 0);
          validNorm = true;
        }
        else
        {
          n = UVector3(fSinSPhi, -fCosSPhi, 0);
          validNorm = false;
        }
        break;

      case kEPhi:
        if (fDPhi <= UUtils::kPi)
        {
          n = UVector3(-fSinEPhi, fCosEPhi, 0);
          validNorm = true;
        }
        else
        {
          n = UVector3(-fSinEPhi, fCosEPhi, 0);
          validNorm = false;
        }
        break;

      case kPZ:
        n        = UVector3(0, 0, 1);
        validNorm = true;
        break;

      case kMZ:
        n        = UVector3(0, 0, -1);
        validNorm = true;
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
                << "Direction:" << std::endl << std::endl
                << "v.x = "  << v.x() << std::endl
                << "v.y = "  << v.y() << std::endl
                << "v.z = "  << v.z() << std::endl << std::endl
                << "Proposed distance :" << std::endl << std::endl
                << "snxt = "    << snxt << " mm" << std::endl;
        message.precision(oldprc);
        UUtils::Exception("UTubs::DistanceToOut(p,v,..)", "GeomSolids1002",
                          UWarning, 1, message.str().c_str());
        break;
    }
  }
  if (snxt < halfCarTolerance)
  {
    snxt = 0;
  }

  return snxt;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

double UTubs::SafetyFromInside(const UVector3& p, bool) const
{
  double safe = 0.0, rho, safeZ;
  rho = std::sqrt(p.x() * p.x() + p.y() * p.y());

#ifdef UDEBUG
  if (Inside(p) == eOutside)
  {
    int oldprc = cout.precision(16);
    cout << std::endl;
    DumpInfo();
    cout << "Position:" << std::endl << std::endl;
    cout << "p.x = "   << p.x() << " mm" << std::endl;
    cout << "p.y = "   << p.y() << " mm" << std::endl;
    cout << "p.z = "   << p.z() << " mm" << std::endl << std::endl;
    cout.precision(oldprc);
    UUtils::Exception("UTubs::DistanceToOut(p)", "GeomSolids1002",
                      UWarning, 1, "Point p is outside !?");
  }
#endif
  safe = SafetyFromInsideR(p, rho);
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
// Stream object contents to an output stream

UGeometryType UTubs::GetEntityType() const
{
  return std::string("Tubs");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
VUSolid* UTubs::Clone() const
{
  return new UTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& UTubs::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "                *** Dump for solid - " << GetName() << " ***\n"
     << "                ===================================================\n"
     << " Solid type: UTubs\n"
     << " Parameters: \n"
     << "                inner radius : " << fRMin << " mm \n"
     << "                outer radius : " << fRMax << " mm \n"
     << "                half length Z: " << fDz << " mm \n"
     << "                starting phi : " << fSPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "                delta phi    : " << fDPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

UVector3 UTubs::GetPointOnSurface() const
{
  double xRand, yRand, zRand, phi, cosphi, sinphi, chose,
         aOne, aTwo, aThr, aFou;
  double rRand;

  aOne = 2.*fDz * fDPhi * fRMax;
  aTwo = 2.*fDz * fDPhi * fRMin;
  aThr = 0.5 * fDPhi * (fRMax * fRMax - fRMin * fRMin);
  aFou = 2.*fDz * (fRMax - fRMin);

  phi   = UUtils::Random(fSPhi, fSPhi + fDPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  rRand = UUtils::GetRadiusInRing(fRMin, fRMax);

  if ((fSPhi == 0) && (fDPhi == 2 * UUtils::kPi))
  {
    aFou = 0;
  }

  chose = UUtils::Random(0., aOne + aTwo + 2.*aThr + 2.*aFou);

  if ((chose >= 0) && (chose < aOne))
  {
    xRand = fRMax * cosphi;
    yRand = fRMax * sinphi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand);
  }
  else if ((chose >= aOne) && (chose < aOne + aTwo))
  {
    xRand = fRMin * cosphi;
    yRand = fRMin * sinphi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand);
  }
  else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aThr))
  {
    xRand = rRand * cosphi;
    yRand = rRand * sinphi;
    zRand = fDz;
    return UVector3(xRand, yRand, zRand);
  }
  else if ((chose >= aOne + aTwo + aThr) && (chose < aOne + aTwo + 2.*aThr))
  {
    xRand = rRand * cosphi;
    yRand = rRand * sinphi;
    zRand = -1.*fDz;
    return UVector3(xRand, yRand, zRand);
  }
  else if ((chose >= aOne + aTwo + 2.*aThr)
           && (chose < aOne + aTwo + 2.*aThr + aFou))
  {
    xRand = rRand * fCosSPhi;
    yRand = rRand * fSinSPhi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand);
  }
  else
  {
    xRand = rRand * fCosSPhiDPhi;
    yRand = rRand * fSinSPhiDPhi;
    zRand = UUtils::Random(-1.*fDz, fDz);
    return UVector3(xRand, yRand, zRand);
  }
}

void UTubs::Extent(UVector3& aMin, UVector3& aMax) const
{
  aMin = UVector3(-fRMax, -fRMax, -fDz);
  aMax = UVector3(fRMax, fRMax, fDz);
}

void UTubs::GetParametersList(int, double* aArray) const
{
  aArray[0] = GetInnerRadius();
  aArray[1] = GetOuterRadius();
  aArray[2] = GetZHalfLength();
  aArray[3] = GetStartPhiAngle();
  aArray[4] = GetDeltaPhiAngle();
}
