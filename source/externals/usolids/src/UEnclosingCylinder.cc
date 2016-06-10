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
// UEnclosingCylinder
//
// 19.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UEnclosingCylinder.hh"
#include "UReduciblePolygon.hh"

#include "VUSolid.hh"
#include "UBox.hh"
#include "UTubs.hh"

//
// Constructor
//
UEnclosingCylinder::UEnclosingCylinder(/*const UReduciblePolygon *rz*/double r, double hi, double lo,
    bool thePhiIsOpen,
    double theStartPhi,
    double theTotalPhi)
  : startPhi(theStartPhi), totalPhi(theTotalPhi),
    rx1(0.), ry1(0.), dx1(0.), dy1(0.),
    rx2(0.), ry2(0.), dx2(0.), dy2(0.),
    concave(theTotalPhi > UUtils::kPi)
{
  //
  // Obtain largest r and smallest and largest z
  //

  /*
  radius = rz->Amax();
  zHi = rz->Bmax();
  zLo = rz->Bmin();
  */

  radius = r;
  zHi = hi;
  zLo = lo;

  double fTolerance = VUSolid::Tolerance();

  tube = new UTubs("", 0, radius + fTolerance, zHi - zLo, theStartPhi, theTotalPhi);

  //
  // Save phi info
  //
  phiIsOpen = thePhiIsOpen;
  if (phiIsOpen)
  {
    rx1 = std::cos(startPhi);
    ry1 = std::sin(startPhi);
    dx1 = +ry1 * 10 * fTolerance;
    dy1 = -rx1 * 10 * fTolerance;

    rx2 = std::cos(startPhi + totalPhi);
    ry2 = std::sin(startPhi + totalPhi);
    dx2 = -ry2 * 10 * fTolerance;
    dy2 = +rx2 * 10 * fTolerance;
  }

  //
  // Add safety
  //
  radius += 10 * fTolerance;
  zLo   -= 10 * fTolerance;
  zHi   += 10 * fTolerance;
}


//
// Destructor
//
UEnclosingCylinder::~UEnclosingCylinder()
{
}


//
// Outside
//
// Decide very rapidly if the point is outside the cylinder
//
// If one is not certain, return false
//
bool UEnclosingCylinder::MustBeOutside(const UVector3& p) const
{
//  if (p.Perp() > radius ) return true;
  if (p.Perp2() > radius * radius) return true;
  if (p.z() < zLo) return true;
  if (p.z() > zHi) return true;

  if (phiIsOpen)
  {
    if (concave)
    {
      if (((p.x() - dx1)*ry1 - (p.y() - dy1)*rx1) < 0) return false;
      if (((p.x() - dx2)*ry2 - (p.y() - dy2)*rx2) > 0) return false;
    }
    else
    {
      if (((p.x() - dx1)*ry1 - (p.y() - dy1)*rx1) > 0) return true;
      if (((p.x() - dx2)*ry2 - (p.y() - dy2)*rx2) < 0) return true;
    }
  }

  return false;
}


//
// Misses
//
// Decide very rapidly if the trajectory is going to miss the cylinder
//
// If one is not sure, return false
//
bool UEnclosingCylinder::ShouldMiss(const UVector3& p,
                                    const UVector3& v) const
{
  if (!MustBeOutside(p)) return false;

//  if (p.z() < zLo - VUSolid::Tolerance() && v.z() <= 0) return true;
//  if (p.z() > zHi + VUSolid::Tolerance() && v.z() >= 0) return true;

  double cross = p.x() * v.y() - p.y() * v.x();
  if (cross > radius) return true;

  double r2 = p.Perp2();
  if (r2 > radius * radius)
  {
    double dot = p.x() * v.x() + p.y() * v.y();
    if (dot > 0) return true;

//    double n = v.Perp2();
//    if (Dot < std::sqrt(r2 - radius * radius) * n) return true;
  }


  /*
  if (phiIsOpen)
  {
    if (concave)
    {
      if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) < 0) return false;
      if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) > 0) return false;
    }
    else
    {
      if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) > 0) return true;
      if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) < 0) return true;
    }
    return false;
  }
  */

  return false;
}

void UEnclosingCylinder::Extent(UVector3& aMin, UVector3& aMax) const
{
  aMin = UVector3(-radius, -radius, zLo);
  aMax = UVector3(radius, radius, zHi);
}

double UEnclosingCylinder::DistanceTo(const UVector3& p, const UVector3& v) const
{
  return tube->DistanceToIn(p, v);
}

double UEnclosingCylinder::SafetyFromOutside(const UVector3& p) const
{
  return tube->SafetyFromOutside(p);
}
