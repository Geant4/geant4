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
// UQuadrangularFacet
//
// 17.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include "VUSolid.hh"
#include "UQuadrangularFacet.hh"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//
// !!!THIS IS A FUDGE!!!  IT'S TWO ADJACENT UTRIANGULARFACETS
// --- NOT EFFICIENT BUT PRACTICAL.
//
UQuadrangularFacet::UQuadrangularFacet(const UVector3& vt0,
                                       const UVector3& vt1, const UVector3& vt2,
                                       const UVector3& vt3, UFacetVertexType vertexType)
{
  UVector3 e1, e2, e3;

  SetVertex(0, vt0);
  if (vertexType == UABSOLUTE)
  {
    SetVertex(1, vt1);
    SetVertex(2, vt2);
    SetVertex(3, vt3);

    e1 = vt1 - vt0;
    e2 = vt2 - vt0;
    e3 = vt3 - vt0;
  }
  else
  {
    SetVertex(1, vt0 + vt1);
    SetVertex(2, vt0 + vt2);
    SetVertex(3, vt0 + vt3);

    e1 = vt1;
    e2 = vt2;
    e3 = vt3;
  }
  double length1 = e1.Mag();
  double length2 = (GetVertex(2) - GetVertex(1)).Mag();
  double length3 = (GetVertex(3) - GetVertex(2)).Mag();
  double length4 = e3.Mag();

  UVector3 normal1 = e1.Cross(e2).Unit();
  UVector3 normal2 = e2.Cross(e3).Unit();

  bool isDefined = (length1 > VUSolid::Tolerance() && length2 > VUSolid::Tolerance() &&
                    length3 > VUSolid::Tolerance() && length4 > VUSolid::Tolerance() &&
                    normal1.Dot(normal2) >= 0.9999999999);

  if (isDefined)
  {
    fFacet1 = UTriangularFacet(GetVertex(0), GetVertex(1), GetVertex(2), UABSOLUTE);
    fFacet2 = UTriangularFacet(GetVertex(0), GetVertex(2), GetVertex(3), UABSOLUTE);

    UTriangularFacet facet3(GetVertex(0), GetVertex(1), GetVertex(3), UABSOLUTE);
    UTriangularFacet facet4(GetVertex(1), GetVertex(2), GetVertex(3), UABSOLUTE);

    UVector3 normal12 = fFacet1.GetSurfaceNormal() + fFacet2.GetSurfaceNormal();
    UVector3 normal34 = facet3.GetSurfaceNormal() + facet4.GetSurfaceNormal();
    UVector3 normal = 0.25 * (normal12 + normal34);

    fFacet1.SetSurfaceNormal(normal);
    fFacet2.SetSurfaceNormal(normal);

    UVector3 vtmp = 0.5 * (e1 + e2);
    fCircumcentre = GetVertex(0) + vtmp;
    double radiusSqr = vtmp.Mag2();
    fRadius = std::sqrt(radiusSqr);
  }
  else
  {
    UUtils::Exception("UQuadrangularFacet::UQuadrangularFacet()", "GeomSolids1002",
        UWarning, 1, "Length of sides of facet are too small or sides not planar.");
    cerr << endl;
    cerr << "P0 = " << GetVertex(0) << endl;
    cerr << "P1 = " << GetVertex(1) << endl;
    cerr << "P2 = " << GetVertex(2) << endl;
    cerr << "P3 = " << GetVertex(3) << endl;
    cerr << "Side lengths = P0->P1" << length1 << endl;
    cerr << "Side lengths = P1->P2" << length2 << endl;
    cerr << "Side lengths = P2->P3" << length3 << endl;
    cerr << "Side lengths = P3->P0" << length4 << endl;
    cerr << endl;
    fRadius = 0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
UQuadrangularFacet::~UQuadrangularFacet()
{
}

///////////////////////////////////////////////////////////////////////////////
//
UQuadrangularFacet::UQuadrangularFacet(const UQuadrangularFacet& rhs) : VUFacet(rhs)
{
  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;
  fRadius = 0;
}

///////////////////////////////////////////////////////////////////////////////
//
UQuadrangularFacet& UQuadrangularFacet::operator=(const UQuadrangularFacet& rhs)
{
  if (this == &rhs)
    return *this;

  fFacet1 = rhs.fFacet1;
  fFacet2 = rhs.fFacet2;

  fRadius = 0;

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
VUFacet* UQuadrangularFacet::GetClone()
{
  UQuadrangularFacet* c = new UQuadrangularFacet(GetVertex(0), GetVertex(1), GetVertex(2), GetVertex(3), UABSOLUTE);
  return c;
}

///////////////////////////////////////////////////////////////////////////////
//
UVector3 UQuadrangularFacet::Distance(const UVector3& p)
{
  UVector3 v1 = fFacet1.Distance(p);
  UVector3 v2 = fFacet2.Distance(p);

  if (v1.Mag2() < v2.Mag2()) return v1;
  else return v2;
}

///////////////////////////////////////////////////////////////////////////////
//
double UQuadrangularFacet::Distance(const UVector3& p,
                                    const double)
{
  double dist = Distance(p).Mag();
  return dist;
}

///////////////////////////////////////////////////////////////////////////////
//
double UQuadrangularFacet::Distance(const UVector3& p, const double, const bool outgoing)
{
  double dist;

  UVector3 v = Distance(p);
  double dir = v.Dot(GetSurfaceNormal());
  if ((dir > dirTolerance && !outgoing) || (dir < -dirTolerance && outgoing))
    dist = UUtils::kInfinity;
  else
    dist = v.Mag();
  return dist;
}

double UQuadrangularFacet::Extent(const UVector3 axis)
{
  double ss  = 0;

  for (int i = 0; i <= 3; ++i)
  {
    double sp = GetVertex(i).Dot(axis);
    if (sp > ss) ss = sp;
  }
  return ss;
}

bool UQuadrangularFacet::Intersect(const UVector3& p, const UVector3& v, bool outgoing, double& distance, double& distFromSurface, UVector3& normal)
{
  bool intersect = fFacet1.Intersect(p, v, outgoing, distance, distFromSurface, normal);
  if (!intersect) intersect = fFacet2.Intersect(p, v, outgoing, distance, distFromSurface, normal);
  if (!intersect)
  {
    distance = distFromSurface = UUtils::kInfinity;
    normal.Set(0);
  }
  return intersect;
}

// Auxiliary method for get a random point on surface
UVector3 UQuadrangularFacet::GetPointOnFace() const
{
  UVector3 pr = (UUtils::Random(0., 1.) < 0.5) ? fFacet1.GetPointOnFace() : fFacet2.GetPointOnFace();
  return pr;
}

// Auxiliary method for returning the surface area
double UQuadrangularFacet::GetArea()
{
  double area = fFacet1.GetArea() + fFacet2.GetArea();
  return area;
}

std::string UQuadrangularFacet::GetEntityType() const
{
  return "QuadrangularFacet";
}

UVector3 UQuadrangularFacet::GetSurfaceNormal() const
{
  return fFacet1.GetSurfaceNormal();
}
