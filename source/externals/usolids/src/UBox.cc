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
// UBox
//
// 10.06.11 J.Apostolakis, G.Cosmo, A.Gheata
//          Created from original implementation in Geant4 and ROOT
// --------------------------------------------------------------------

#include <iostream>
#include <sstream>
#include "UUtils.hh"
#include "UBox.hh"

//______________________________________________________________________________
UBox::UBox(const std::string& name, double dx, double dy, double dz)
  : VUSolid(name),
    fDx(dx),
    fDy(dy),
    fDz(dz), fCubicVolume(0.), fSurfaceArea(0.)
{
  // Named constructor

  if ((dx < 2 * VUSolid::fgTolerance)
      || (dy < 2 * VUSolid::fgTolerance)
      || (dz < 2 * VUSolid::fgTolerance))  // limit to thickness of surfaces
  {
    //std::ostringstream message;
    std::ostringstream message;
    message << "Dimensions too small for Solid: " << GetName() << "!" << std::endl
            << "     dx, dy, dz = " << dx << ", " << dy << ", " << dz;
    UUtils::Exception("UBox::UBox()", "GeomSolids0002", UFatalErrorInArguments, 1, message.str().c_str());
  }
}

void UBox::Set(double dx, double dy, double dz)
{
  fDx = dx;
  fDy = dy;
  fDz = dz;
}

void UBox::Set(const UVector3& vec)
{
  fDx = vec.x();
  fDy = vec.y();
  fDz = vec.z();
}
//Destructor
UBox::~UBox()
{

}
// Copy constructor

UBox::UBox(const UBox& rhs)
  : VUSolid(rhs), fDx(rhs.fDx), fDy(rhs.fDy), fDz(rhs.fDz),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea)
{
}

//
// Assignment operator

UBox& UBox::operator = (const UBox& rhs)
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
  fDx = rhs.fDx;
  fDy = rhs.fDy;
  fDz = rhs.fDz;
  fCubicVolume = rhs.fCubicVolume;
  fSurfaceArea = rhs.fSurfaceArea;

  return *this;
}

//______________________________________________________________________________
VUSolid::EnumInside UBox::Inside(const UVector3& aPoint) const
{
  // Classify point location with respect to solid:
  //  o eInside       - inside the solid
  //  o eSurface      - close to surface within tolerance
  //  o eOutside      - outside the solid
  static const double delta = VUSolid::fgTolerance;
  // Early returns on outside condition on any axis. Check Z first for faster
  // exclusion in  phi symmetric geometries.
  double ddz = std::abs(aPoint.z()) - fDz;
  if (ddz > delta) return eOutside;
  double ddx = std::abs(aPoint.x()) - fDx;
  if (ddx > delta) return eOutside;
  double ddy = std::abs(aPoint.y()) - fDy;
  if (ddy > delta) return eOutside;
  if (ddx > - delta || ddy > -delta || ddz > -delta) return eSurface;
  return eInside;
}

//______________________________________________________________________________
double UBox::DistanceToIn(const UVector3& aPoint,
                          const UVector3& aDirection,
                          //                          UVector3 &aNormal,
                          double aPstep) const
{
  // Computes distance from a point presumably outside the solid to the solid
  // surface. Ignores first surface if the point is actually inside. Early return
  // infinity in case the safety to any surface is found greater than the proposed
  // step aPstep.
  // The normal vector to the crossed surface is filled only in case the box is
  // crossed, otherwise aNormal.IsNull() is true.

  // Compute safety to the closest surface on each axis.
  // Early exits if safety bigger than proposed step.
  static const double delta = VUSolid::fgTolerance;
  //   aNormal.SetNull();
  double safx = std::abs(aPoint.x()) - fDx;
  double safy = std::abs(aPoint.y()) - fDy;
  double safz = std::abs(aPoint.z()) - fDz;
  if ((safx > aPstep) || (safy > aPstep) || (safz > aPstep))
    return UUtils::kInfinity;
  // Check numerical outside.
  bool outside = (safx > 0) || (safy > 0) || (safz > 0);
  if (!outside)
  {
    // If point close to this surface, check against the normal
    if (safx > -delta)
    {
      if (aPoint.x() * aDirection.x() > 0) return UUtils::kInfinity ;
    }
    if (safy > -delta)
    {
      if (aPoint.y() * aDirection.y() > 0) return UUtils::kInfinity ;
    }
    if (safz > -delta)
    {
      if (aPoint.z() * aDirection.z() > 0) return UUtils::kInfinity ;
    }
    // Point actually "deep" inside, return zero distance, normal un-defined
    return 0.0;
  }
  // The point is really outside. Only axis with positive safety to be
  // considered. Early exit on each axis if point and direction components
  // have the same .sign.
  double dist = 0.0;
  double coordinate = 0.0;
  if (safx > 0)
  {
    if (aPoint.x() * aDirection.x() >= 0) return UUtils::kInfinity;
    dist = safx / std::abs(aDirection.x());
    coordinate = aPoint.y() + dist * aDirection.y();
    if (std::abs(coordinate) < fDy)
    {
      coordinate = aPoint.z() + dist * aDirection.z();
      if (std::abs(coordinate) < fDz)
      {
        //            aNormal.x() = UUtils::Sign(1.0, aPoint.x());
        if (dist < 0.5 * delta) dist = 0.;
        return dist;
      }
    }
  }
  if (safy > 0)
  {
    if (aPoint.y() * aDirection.y() >= 0) return UUtils::kInfinity;
    dist = safy / std::abs(aDirection.y());
    coordinate = aPoint.x() + dist * aDirection.x();
    if (std::abs(coordinate) < fDx)
    {
      coordinate = aPoint.z() + dist * aDirection.z();
      if (std::abs(coordinate) < fDz)
      {
        //            aNormal.y() = UUtils::Sign(1.0, aPoint.y());
        if (dist < 0.5 * delta) dist = 0.;
        return dist;
      }
    }
  }
  if (safz > 0)
  {
    if (aPoint.z() * aDirection.z() >= 0) return UUtils::kInfinity;
    dist = safz / std::abs(aDirection.z());
    coordinate = aPoint.x() + dist * aDirection.x();
    if (std::abs(coordinate) < fDx)
    {
      coordinate = aPoint.y() + dist * aDirection.y();
      if (std::abs(coordinate) < fDy)
      {
        //            aNormal.z() = UUtils::Sign(1.0, aPoint.z());
        if (dist < 0.5 * delta) dist = 0.;
        return dist;
      }
    }
  }
  return UUtils::kInfinity;
}

//______________________________________________________________________________
double UBox::DistanceToOut(const UVector3&  aPoint, const UVector3& aDirection,
                           UVector3& aNormal,
                           bool&    convex,
                           double /*aPstep*/) const
{
  // Computes distance from a point presumably intside the solid to the solid
  // surface. Ignores first surface along each axis systematically (for points
  // inside or outside. Early returns zero in case the second surface is behind
  // the starting point.
  // o The proposed step is ignored.
  // o The normal vector to the crossed surface is always filled.
  double smin = UUtils::kInfinity;
  double snxt, signDir;
  convex = true;  //   Box is convex (even if the starting point is outside)
  // Check always the "away" surface along direction on axis. This responds
  // corectly even for points outside the solid (no need for tolerance check)
  if (aDirection.x() != 0.0)
    { 
    signDir = UUtils::Sign(1.0, aDirection.x());
    aNormal.Set(signDir, 0., 0.);
    snxt = (-aPoint.x() + signDir * fDx) / aDirection.x();
    if (snxt <= 0) return 0.0;   // point outside moving outwards
    smin = snxt;
  }

  if (aDirection.y() != 0.0)
    { 
    signDir = UUtils::Sign(1.0, aDirection.y());
    snxt = (-aPoint.y() + signDir * fDy) / aDirection.y();
    if (snxt <= 0)
    {
      aNormal.Set(0., signDir, 0.);
      return 0.0; // point outside moving outwards
    }
    if (snxt < smin)
    {
      smin = snxt;
      aNormal.Set(0., signDir, 0.);
    }
  }

  if (aDirection.z() != 0.0)
    { 
    signDir = UUtils::Sign(1.0, aDirection.z());
    snxt = (-aPoint.z() + signDir * fDz) / aDirection.z();
    if (snxt <= 0)
    {
      aNormal.Set(0., 0., signDir);
      return 0.0; // point outside moving outwards
    }
    if (snxt < smin)
    {
      smin = snxt;
      aNormal.Set(0., 0., signDir);
     }
  }
  if (smin < 0.5 * VUSolid::fgTolerance) smin = 0.;
  return smin;
}

//______________________________________________________________________________
double UBox::SafetyFromInside(const UVector3& aPoint,
                              bool /*aAccurate*/) const
{
  // Estimates the isotropic safety from a point inside the current solid to any
  // of its surfaces. The algorithm may be accurate or should provide a fast
  // underestimate.
  double safe, safy, safz;
  safe = fDx - std::abs(aPoint.x());
  safy = fDy - std::abs(aPoint.y());
  if (safy < safe) safe = safy;
  safz = fDz - std::abs(aPoint.z());
  if (safz < safe) safe = safz;
  return std::max(0.0, safe);
}

//______________________________________________________________________________
double UBox::SafetyFromOutside(const UVector3& aPoint,
                               bool aAccurate) const
{
  // Estimates the isotropic safety from a point outside the current solid to any
  // of its surfaces. The algorithm may be accurate or should provide a fast
  // underestimate.
  double safe, safx, safy, safz;
  safe = safx = -fDx + std::abs(aPoint.x());
  safy = -fDy + std::abs(aPoint.y());
  if (safy > safe) safe = safy;
  safz = -fDz + std::abs(aPoint.z());
  if (safz > safe) safe = safz;
  if (safe < 0.0) return 0.0; // point is inside
  if (!aAccurate) return safe;
  double safsq = 0.0;
  int count = 0;
  if (safx > 0)
  {
    safsq += safx * safx;
    count++;
  }
  if (safy > 0)
  {
    safsq += safy * safy;
    count++;
  }
  if (safz > 0)
  {
    safsq += safz * safz;
    count++;
  }
  if (count == 1) return safe;
  return std::sqrt(safsq);
}

//______________________________________________________________________________
bool UBox::Normal(const UVector3& aPoint, UVector3& aNormal) const
{
  // Computes the normal on a surface and returns it as a unit vector
  //   In case a point is further than tolerance_normal from a surface, set validNormal=false
  //   Must return a valid vector. (even if the point is not on the surface.)
  //
  //   On an edge or corner, provide an average normal of all facets within tolerance
  // NOTE: the tolerance value used in here is not yet the global surface
  //     tolerance - we will have to revise this value - TODO
  static const double delta = 100.*VUSolid::fgTolerance;
  static const double kInvSqrt2 = 1. / std::sqrt(2.);
  static const double kInvSqrt3 = 1. / std::sqrt(3.);
  aNormal.Set(0.);
  UVector3 crt_normal, min_normal;
  int nsurf = 0;
  double safx = std::abs(std::abs(aPoint.x()) - fDx);
  double safmin = safx;
  crt_normal.Set(UUtils::Sign(1., aPoint.x()), 0., 0.);
  min_normal = crt_normal;
  if (safx < delta)
  {
    nsurf++;
    aNormal += crt_normal;
  }
  double safy = std::abs(std::abs(aPoint.y()) - fDy);
  crt_normal.Set(0., UUtils::Sign(1., aPoint.y()), 0.);
  if (safy < delta)
  {
    nsurf++;
    aNormal += crt_normal;
  }
  if (safy < safmin)
  {
    min_normal = crt_normal;
    safmin = safy;
  }
  double safz = std::abs(std::abs(aPoint.z()) - fDz);
  crt_normal.Set(0., 0., UUtils::Sign(1., aPoint.z()));
  if (safz < delta)
  {
    nsurf++;
    aNormal += crt_normal;
  }
  if (safz < safmin)
  {
    min_normal = crt_normal;
    safmin = safz;
  }

  bool valid = true;
  switch (nsurf)
  {
    case 0:
      aNormal = min_normal;
      valid = false;
      break;
    case 1:
      break;
    case 2:
      aNormal *= kInvSqrt2;
      break;
    case 3:
      aNormal *= kInvSqrt3;
  };
  return valid;
}


//______________________________________________________________________________
void UBox::Extent(UVector3& aMin, UVector3& aMax) const
{
  // Returns the full 3D cartesian extent of the solid.
  aMin.x() = -fDx;
  aMax.x() = fDx;
  aMin.y() = -fDy;
  aMax.y() = fDy;
  aMin.z() = -fDz;
  aMax.z() = fDz;
}

/////////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Return a point (UVector3) randomly and uniformly selected
// on the solid surface

UVector3 UBox::GetPointOnSurface() const
{
  double px, py, pz, select, sumS;
  double Sxy = fDx * fDy, Sxz = fDx * fDz, Syz = fDy * fDz;

  sumS   = Sxy + Sxz + Syz;
  select = sumS * UUtils::Random();

  if (select < Sxy)
  {
    px = -fDx + 2 * fDx * UUtils::Random();
    py = -fDy + 2 * fDy * UUtils::Random();

    if (UUtils::Random() > 0.5)
    {
      pz = fDz;
    }
    else
    {
      pz = -fDz;
    }
  }
  else if ((select - Sxy) < Sxz)
  {
    px = -fDx + 2 * fDx * UUtils::Random();
    pz = -fDz + 2 * fDz * UUtils::Random();

    if (UUtils::Random() > 0.5)
    {
      py = fDy;
    }
    else
    {
      py = -fDy;
    }
  }
  else
  {
    py = -fDy + 2 * fDy * UUtils::Random();
    pz = -fDz + 2 * fDz * UUtils::Random();

    if (UUtils::Random() > 0.5)
    {
      px = fDx;
    }
    else
    {
      px = -fDx;
    }
  }
  return UVector3(px, py, pz);
}

/////////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnEdge
//
// Return a point (UVector3) randomly and uniformly selected
// on the solid edge

UVector3 UBox::GetPointOnEdge() const
{
  double select, sumL;
  double Lx = 2 * fDx, Ly = 2 * fDy, Lz= 2 * fDz;

  sumL   = Lx + Ly + Lz;
  select = sumL * UUtils::Random();

  if (select < Lx)
  {
    select = UUtils::Random();
    if (select < 0.25) return UVector3( -fDx + 2 * fDx* UUtils::Random(),-fDy,-fDz);
    if (select < 0.5 ) return UVector3( -fDx + 2 * fDx* UUtils::Random(), fDy,-fDz);
    if (select < 0.75) return UVector3( -fDx + 2 * fDx* UUtils::Random(),-fDy, fDz);
    return UVector3( -fDx + 2 * fDx* UUtils::Random(),fDy,fDz);
  }
  else if ((select - Lx) < Ly)
  {
    select = UUtils::Random();
    if (select < 0.25) return UVector3( -fDx,-fDy + 2 * fDy* UUtils::Random(),-fDz);
    if (select < 0.5 ) return UVector3(  fDx,-fDy + 2 * fDy* UUtils::Random(),-fDz);
    if (select < 0.75) return UVector3( -fDx,-fDy + 2 * fDy* UUtils::Random(), fDz);
    return UVector3( fDx,-fDy + 2 * fDy* UUtils::Random(),fDz);
  }
  else
  {
   select = UUtils::Random();
   if (select < 0.25) return UVector3( -fDx,-fDy,-fDz + 2 * fDz* UUtils::Random());
   if (select < 0.5 ) return UVector3(  fDx,-fDy,-fDz + 2 * fDz* UUtils::Random());
   if (select < 0.75) return UVector3( -fDx, fDy,-fDz + 2 * fDz* UUtils::Random());
   return UVector3( fDx,fDy,-fDz + 2 * fDz* UUtils::Random());
    
  }
  return 0;
}

std::ostream& UBox::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "		*** Dump for solid - " << GetName() << " ***\n"
     << "		===================================================\n"
     << " Solid type: UBox\n"
     << " Parameters: \n"
     << "		half length X: " << fDx << " mm \n"
     << "		half length Y: " << fDy << " mm \n"
     << "		half length Z: " << fDz << " mm \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

void UBox::SetXHalfLength(double dx)
{
  if (dx > 2 * VUSolid::fgTolerance) // limit to thickness of surfaces
  {
    fDx = dx;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension X too small for solid: " << GetName() << "!"
            << std::endl
            << "       hX = " << dx;
    UUtils::Exception("UBox::SetXHalfLength()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;

}

void UBox::SetYHalfLength(double dy)
{
  if (dy > 2 * VUSolid::fgTolerance) // limit to thickness of surfaces
  {
    fDy = dy;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Y too small for solid: " << GetName() << "!"
            << std::endl
            << "       hY = " << dy;
    UUtils::Exception("UBox::SetYHalfLength()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;

}

void UBox::SetZHalfLength(double dz)
{
  if (dz > 2 * VUSolid::fgTolerance) // limit to thickness of surfaces
  {
    fDz = dz;
  }
  else
  {
    std::ostringstream message;
    message << "Dimension Z too small for solid: " << GetName() << "!"
            << std::endl
            << "       hZ = " << dz;
    UUtils::Exception("G4Box::SetZHalfLength()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;

}

UGeometryType UBox::GetEntityType() const
{
   return "Box";
}

VUSolid* UBox::Clone() const
{
   return new UBox(*this);
}
