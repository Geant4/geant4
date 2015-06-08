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
// UExtrudedSolid
//
// 13.08.13 Tatiana Nikitina
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include <iomanip>
#include <sstream>
#include <set>

#include "UExtrudedSolid.hh"
#include "VUFacet.hh"
#include "UTriangularFacet.hh"
#include "UQuadrangularFacet.hh"

//_____________________________________________________________________________

UExtrudedSolid::UExtrudedSolid(const std::string& pName,
                               std::vector<UVector2> polygon,
                               std::vector<ZSection> zsections)
  : UTessellatedSolid(pName),
    fPolygon(),
    fZSections(),
    fTriangles(),
    fIsConvex(false),
    fGeometryType("ExtrudedSolid")

{
  // General constructor

  Initialise(polygon, zsections);
}

//_____________________________________________________________________________

UExtrudedSolid::UExtrudedSolid(const std::string& pName,
                               std::vector<UVector2> polygon, double dz,
                               UVector2 off1, double scale1,
                               UVector2 off2, double scale2)
  : UTessellatedSolid(pName),
    fNz(2),
    fPolygon(),
    fZSections(),
    fTriangles(),
    fIsConvex(false),
    fGeometryType("ExtrudedSolid")

{
  // Special constructor for solid with 2 z-sections

  Initialise(polygon, dz, off1, scale1, off2, scale2);
}

//_____________________________________________________________________________

UExtrudedSolid::UExtrudedSolid()
  : UTessellatedSolid(), fNv(0), fNz(0), fPolygon(), fZSections(),
    fTriangles(), fIsConvex(false), fGeometryType("UExtrudedSolid")
{
  // Fake default constructor - sets only member data and allocates memory
  //                            for usage restricted to object persistency.
}

//_____________________________________________________________________________

UExtrudedSolid::UExtrudedSolid(const UExtrudedSolid& rhs)
  : UTessellatedSolid(rhs), fNv(rhs.fNv), fNz(rhs.fNz),
    fPolygon(rhs.fPolygon), fZSections(rhs.fZSections),
    fTriangles(rhs.fTriangles), fIsConvex(rhs.fIsConvex),
    fGeometryType(rhs.fGeometryType), fKScales(rhs.fKScales),
    fScale0s(rhs.fScale0s), fKOffsets(rhs.fKOffsets), fOffset0s(rhs.fOffset0s)
{
}


//_____________________________________________________________________________

UExtrudedSolid& UExtrudedSolid::operator = (const UExtrudedSolid& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  UTessellatedSolid::operator=(rhs);

  // Copy data
  //
  fNv = rhs.fNv;
  fNz = rhs.fNz;
  fPolygon = rhs.fPolygon;
  fZSections = rhs.fZSections;
  fTriangles = rhs.fTriangles;
  fIsConvex = rhs.fIsConvex;
  fGeometryType = rhs.fGeometryType;
  fKScales = rhs.fKScales;
  fScale0s = rhs.fScale0s;
  fKOffsets = rhs.fKOffsets;
  fOffset0s = rhs.fOffset0s;

  return *this;
}

//_____________________________________________________________________________

UExtrudedSolid::~UExtrudedSolid()
{
  // Destructor
}

//_____________________________________________________________________________

void UExtrudedSolid::Initialise(std::vector<UVector2>& polygon,
                                std::vector<ZSection>& zsections)
{
  fNv = polygon.size();
  fNz = zsections.size();

  // First check input parameters
  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of polygon vertices < 3 - " << GetName().c_str();
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0002",
                      UFatalErrorInArguments, 2, message.str().c_str());
  }

  if (fNz < 2)
  {
    std::ostringstream message;
    message << "Number of z-sides < 2 - " << GetName().c_str();
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0002",
                      UFatalErrorInArguments, 2, message.str().c_str());
  }

  for (int i = 0; i < fNz - 1; ++i)
  {
    if (zsections[i].fZ > zsections[i + 1].fZ)
    {
      std::ostringstream message;
      message << "Z-sections have to be ordered by z value (z0 < z1 < z2...) - "
              << GetName().c_str();
      UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0002",
                        UFatalErrorInArguments, 2, message.str().c_str());
    }
    if (std::fabs(zsections[i + 1].fZ - zsections[i].fZ) < VUSolid::fgTolerance * 0.5)
    {
      std::ostringstream message;
      message << "Z-sections with the same z position are not supported - "
              << GetName().c_str();
      UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0001",
                        UFatalError, 1, message.str().c_str());
    }
  }

  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)
  //
  double area = 0.;
  for (int i = 0; i < fNv; ++i)
  {
    int j = i + 1;
    if (j == fNv) j = 0;
    area += 0.5 * (polygon[i].x * polygon[j].y - polygon[j].x * polygon[i].y);
  }

  // Copy polygon
  //
  if (area < 0.)
  {
    // Polygon vertices are defined clockwise, we just copy the polygon
    for (int i = 0; i < fNv; ++i)
    {
      fPolygon.push_back(polygon[i]);
    }
  }
  else
  {
    // Polygon vertices are defined anti-clockwise, we revert them
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids1001",
                      UWarning, 4,
                "Polygon vertices defined anti-clockwise, reverting polygon");
    for (int i = 0; i < fNv; ++i)
    {
      fPolygon.push_back(polygon[fNv - i - 1]);
    }
  }


  // Copy z-sections
  //
  for (int i = 0; i < fNz; ++i)
  {
    fZSections.push_back(zsections[i]);
  }


  bool result = MakeFacets();
  if (!result)
  {
    std::ostringstream message;
    message << "Making facets failed - " << GetName().c_str();
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0003",
                     UFatalError, 1, message.str().c_str());
  }
  fIsConvex = IsConvex();

  ComputeProjectionParameters();
}

//_____________________________________________________________________________

void UExtrudedSolid::Initialise(std::vector<UVector2>& polygon, double dz,
                                UVector2 off1, double scale1,
                                UVector2 off2, double scale2)
{

  fNv = polygon.size();

  // First check input parameters
  //
  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of polygon vertices < 3 - " << GetName().c_str();
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0002",
                      UFatalErrorInArguments, 2, message.str().c_str());
  }

  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)

  double area = 0.;
  for (int i = 0; i < fNv; ++i)
  {
    int j = i + 1;
    if (j == fNv)
    {
      j = 0;
    }
    area += 0.5 * (polygon[i].x * polygon[j].y
                   - polygon[j].x * polygon[i].y);
  }

  // Copy polygon
  //
  if (area < 0.)
  {
    // Polygon vertices are defined clockwise, we just copy the polygon
    for (int i = 0; i < fNv; ++i)
    {
      fPolygon.push_back(polygon[i]);
    }
  }
  else
  {
    // Polygon vertices are defined anti-clockwise, we revert them
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids1001",
                      UWarning, 4,
                 "Polygon vertices defined anti-clockwise, reverting polygon");
    for (int i = 0; i < fNv; ++i)
    {
      fPolygon.push_back(polygon[fNv - i - 1]);
    }
  }

  // Copy z-sections
  //
  fZSections.push_back(ZSection(-dz, off1, scale1));
  fZSections.push_back(ZSection(dz, off2, scale2));

  bool result = MakeFacets();
  if (!result)
  {
    std::ostringstream message;
    message << "Making facets failed - " << GetName().c_str();
    UUtils::Exception("UExtrudedSolid::UExtrudedSolid()", "GeomSolids0003",
                      UFatalError, 1, message.str().c_str());
  }
  fIsConvex = IsConvex();

  ComputeProjectionParameters();
}

//_____________________________________________________________________________

void UExtrudedSolid::ComputeProjectionParameters()
{
  // Compute parameters for point projections p(z)
  // to the polygon scale & offset:
  // scale(z) = k*z + scale0
  // offset(z) = l*z + offset0
  // p(z) = scale(z)*p0 + offset(z)
  // p0 = (p(z) - offset(z))/scale(z);
  //

  for (int iz = 0; iz < fNz - 1; ++iz)
  {
    double z1      = fZSections[iz].fZ;
    double z2      = fZSections[iz + 1].fZ;
    double scale1  = fZSections[iz].fScale;
    double scale2  = fZSections[iz + 1].fScale;
    UVector2 off1 = fZSections[iz].fOffset;
    UVector2 off2 = fZSections[iz + 1].fOffset;

    double kscale = (scale2 - scale1) / (z2 - z1);
    double scale0 =  scale2 - kscale * (z2 - z1) / 2.0;
    UVector2 koff = (off2 - off1) / (z2 - z1);
    UVector2 off0 =  off2 - koff * (z2 - z1) / 2.0;

    fKScales.push_back(kscale);
    fScale0s.push_back(scale0);
    fKOffsets.push_back(koff);
    fOffset0s.push_back(off0);
  }
}


//_____________________________________________________________________________

UVector3 UExtrudedSolid::GetVertex(int iz, int ind) const
{
  // Shift and scale vertices

  return UVector3(fPolygon[ind].x * fZSections[iz].fScale
                  + fZSections[iz].fOffset.x,
                  fPolygon[ind].y * fZSections[iz].fScale
                  + fZSections[iz].fOffset.y, fZSections[iz].fZ);
}

//_____________________________________________________________________________


UVector2 UExtrudedSolid::ProjectPoint(const UVector3& point) const
{
  // Project point in the polygon scale
  // scale(z) = k*z + scale0
  // offset(z) = l*z + offset0
  // p(z) = scale(z)*p0 + offset(z)
  // p0 = (p(z) - offset(z))/scale(z);

  // Select projection (z-segment of the solid) according to p.z()
  //
  int iz = 0;
  while (point.z() > fZSections[iz + 1].fZ && iz < fNz - 2)
  {
    ++iz;
  }

  double z0 = (fZSections[iz + 1].fZ + fZSections[iz].fZ) / 2.0;
  UVector2 p2(point.x(), point.y());
  double pscale  = fKScales[iz] * (point.z() - z0) + fScale0s[iz];
  UVector2 poffset = fKOffsets[iz] * (point.z() - z0) + fOffset0s[iz];

  // pscale is always >0 as it is an interpolation between two
  // positive scale values
  //
  return (p2 - poffset) / pscale;
}

//_____________________________________________________________________________

bool UExtrudedSolid::IsSameLine(UVector2 p,
                                UVector2 l1, UVector2 l2) const
{
  // Return true if p is on the line through l1, l2

  if (l1.x == l2.x)
  {
    return std::fabs(p.x - l1.x) < VUSolid::fgTolerance * 0.5;
  }
  double  slope = ((l2.y - l1.y) / (l2.x - l1.x));
  double predy = l1.y +  slope * (p.x - l1.x);
  double dy = p.y - predy;

  // Calculate perpendicular distance
  //
  // double perpD= std::fabs(dy) / std::sqrt( 1 + slope * slope );
  // bool   simpleComp= (perpD<0.5*VUSolid::fgTolerance);

  // Check perpendicular distance vs tolerance 'directly'
  //
  const double tol = 0.5 * VUSolid::fgTolerance ;
  bool    squareComp = (dy * dy < (1 + slope * slope) * tol * tol);

  // return  simpleComp;
  return squareComp;
}

//_____________________________________________________________________________

bool UExtrudedSolid::IsSameLineSegment(UVector2 p,
                                       UVector2 l1, UVector2 l2) const
{
  // Return true if p is on the line through l1, l2 and lies between
  // l1 and l2

  if (p.x < std::min(l1.x, l2.x) - VUSolid::fgTolerance * 0.5 ||
      p.x > std::max(l1.x, l2.x) + VUSolid::fgTolerance * 0.5 ||
      p.y < std::min(l1.y, l2.y) - VUSolid::fgTolerance * 0.5 ||
      p.y > std::max(l1.y, l2.y) + VUSolid::fgTolerance * 0.5)
  {
    return false;
  }

  return IsSameLine(p, l1, l2);
}

//_____________________________________________________________________________

bool UExtrudedSolid::IsSameSide(UVector2 p1, UVector2 p2,
                                UVector2 l1, UVector2 l2) const
{
  // Return true if p1 and p2 are on the same side of the line through l1, l2

  return ((p1.x - l1.x) * (l2.y - l1.y)
          - (l2.x - l1.x) * (p1.y - l1.y))
         * ((p2.x - l1.x) * (l2.y - l1.y)
            - (l2.x - l1.x) * (p2.y - l1.y)) > 0;
}

//_____________________________________________________________________________

bool UExtrudedSolid::IsPointInside(UVector2 a, UVector2 b,
                                   UVector2 c, UVector2 p) const
{
  // Return true if p is inside of triangle abc or on its edges,
  // else returns false

  // Check extent first
  //
  if ((p.x < a.x && p.x < b.x && p.x < c.x) ||
      (p.x > a.x && p.x > b.x && p.x > c.x) ||
      (p.y < a.y && p.y < b.y && p.y < c.y) ||
      (p.y > a.y && p.y > b.y && p.y > c.y)) return false;

  bool inside
    = IsSameSide(p, a, b, c)
      && IsSameSide(p, b, a, c)
      && IsSameSide(p, c, a, b);

  bool onEdge
    = IsSameLineSegment(p, a, b)
      || IsSameLineSegment(p, b, c)
      || IsSameLineSegment(p, c, a);

  return inside || onEdge;
}

//_____________________________________________________________________________

double
UExtrudedSolid::GetAngle(UVector2 po, UVector2 pa, UVector2 pb) const
{
  // Return the angle of the vertex in po

  UVector2 t1 = pa - po;
  UVector2 t2 = pb - po;

  double result = (std::atan2(t1.y, t1.x) - std::atan2(t2.y, t2.x));

  if (result < 0) result += 2 * UUtils::kPi;

  return result;
}

//_____________________________________________________________________________

VUFacet*
UExtrudedSolid::MakeDownFacet(int ind1, int ind2, int ind3) const
{
  // Create a triangular facet from the polygon points given by indices
  // forming the down side ( the normal goes in -z)

  std::vector<UVector3> vertices;
  vertices.push_back(GetVertex(0, ind1));
  vertices.push_back(GetVertex(0, ind2));
  vertices.push_back(GetVertex(0, ind3));

  // first vertex most left
  //
  UVector3 cross
    = (vertices[1] - vertices[0]).Cross(vertices[2] - vertices[1]);

  if (cross.z() > 0.0)
  {
    // vertices ardered clock wise has to be reordered

    UVector3 tmp = vertices[1];
    vertices[1] = vertices[2];
    vertices[2] = tmp;
  }

  return new UTriangularFacet(vertices[0], vertices[1],
                              vertices[2], UABSOLUTE);
}

//_____________________________________________________________________________

VUFacet*
UExtrudedSolid::MakeUpFacet(int ind1, int ind2, int ind3) const
{
  // Creates a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  std::vector<UVector3> vertices;
  vertices.push_back(GetVertex(fNz - 1, ind1));
  vertices.push_back(GetVertex(fNz - 1, ind2));
  vertices.push_back(GetVertex(fNz - 1, ind3));

  // first vertex most left
  //
  UVector3 cross
    = (vertices[1] - vertices[0]).Cross(vertices[2] - vertices[1]);

  if (cross.z() < 0.0)
  {
    // vertices ordered clock wise has to be reordered

    UVector3 tmp = vertices[1];
    vertices[1] = vertices[2];
    vertices[2] = tmp;
  }

  return new UTriangularFacet(vertices[0], vertices[1],
                              vertices[2], UABSOLUTE);
}

//_____________________________________________________________________________

bool UExtrudedSolid::AddGeneralPolygonFacets()
{
  // Decompose polygonal sides in triangular facets

  typedef std::pair < UVector2, int > Vertex;

  // Fill one more vector
  //
  std::vector< Vertex > verticesToBeDone;
  for (int i = 0; i < fNv; ++i)
  {
    verticesToBeDone.push_back(Vertex(fPolygon[i], i));
  }
  std::vector< Vertex > ears;

  std::vector< Vertex >::iterator c1 = verticesToBeDone.begin();
  std::vector< Vertex >::iterator c2 = c1 + 1;
  std::vector< Vertex >::iterator c3 = c1 + 2;
  while (verticesToBeDone.size() > 2)
  {
    // skip concave vertices
    //
    double angle = GetAngle(c2->first, c3->first, c1->first);

    int counter = 0;
    while (angle >= UUtils::kPi)
    {
      // try next three consecutive vertices
      //
      c1 = c2;
      c2 = c3;
      ++c3;
      if (c3 == verticesToBeDone.end())
      {
        c3 = verticesToBeDone.begin();
      }

      angle = GetAngle(c2->first, c3->first, c1->first);

      counter++;

      if (counter > fNv)
      {
        UUtils::Exception("UExtrudedSolid::AddGeneralPolygonFacets",
                          "GeomSolids0003", UFatalError, 1,
                          "Triangularisation has failed.");
        break;
      }
    }

    bool good = true;
    std::vector< Vertex >::iterator it;
    for (it = verticesToBeDone.begin(); it != verticesToBeDone.end(); ++it)
    {
      // skip vertices of tested triangle
      //
      if (it == c1 || it == c2 || it == c3)
      {
        continue;
      }

      if (IsPointInside(c1->first, c2->first, c3->first, it->first))
      {
        good = false;

        // try next three consecutive vertices
        //
        c1 = c2;
        c2 = c3;
        ++c3;
        if (c3 == verticesToBeDone.end())
        {
          c3 = verticesToBeDone.begin();
        }
        break;
      }
    }
    if (good)
    {
      // all points are outside triangle, we can make a facet

      bool result;
      result = AddFacet(MakeDownFacet(c1->second, c2->second, c3->second));
      if (! result)
      {
        return false;
      }

      result = AddFacet(MakeUpFacet(c1->second, c2->second, c3->second));
      if (! result)
      {
        return false;
      }

      std::vector<int> triangle(3);
      triangle[0] = c1->second;
      triangle[1] = c2->second;
      triangle[2] = c3->second;
      fTriangles.push_back(triangle);

      // remove the ear point from verticesToBeDone
      //
      verticesToBeDone.erase(c2);
      c1 = verticesToBeDone.begin();
      c2 = c1 + 1;
      c3 = c1 + 2;
    }
  }
  return true;
}

//_____________________________________________________________________________

bool UExtrudedSolid::MakeFacets()
{
  // Define facets

  bool good;

  // Decomposition of polygonal sides in the facets
  //
  if (fNv == 3)
  {
    good = AddFacet(new UTriangularFacet(GetVertex(0, 0), GetVertex(0, 1),
                                         GetVertex(0, 2), UABSOLUTE));
    if (! good)
    {
      return false;
    }

    good = AddFacet(new UTriangularFacet(GetVertex(fNz - 1, 2), GetVertex(fNz - 1, 1),
                                         GetVertex(fNz - 1, 0), UABSOLUTE));
    if (! good)
    {
      return false;
    }

    std::vector<int> triangle(3);
    triangle[0] = 0;
    triangle[1] = 1;
    triangle[2] = 2;
    fTriangles.push_back(triangle);
  }

  else if (fNv == 4)
  {
    good = AddFacet(new UQuadrangularFacet(GetVertex(0, 0), GetVertex(0, 1),
                                           GetVertex(0, 2), GetVertex(0, 3),
                                           UABSOLUTE));
    if (! good)
    {
      return false;
    }

    good = AddFacet(new UQuadrangularFacet(GetVertex(fNz - 1, 3), GetVertex(fNz - 1, 2),
                                           GetVertex(fNz - 1, 1), GetVertex(fNz - 1, 0),
                                           UABSOLUTE));
    if (! good)
    {
      return false;
    }

    std::vector<int> triangle1(3);
    triangle1[0] = 0;
    triangle1[1] = 1;
    triangle1[2] = 2;
    fTriangles.push_back(triangle1);

    std::vector<int> triangle2(3);
    triangle2[0] = 0;
    triangle2[1] = 2;
    triangle2[2] = 3;
    fTriangles.push_back(triangle2);
  }
  else
  {
    good = AddGeneralPolygonFacets();
    if (! good)
    {
      return false;
    }
  }

  // The quadrangular sides
  //
  for (int iz = 0; iz < fNz - 1; ++iz)
  {
    for (int i = 0; i < fNv; ++i)
    {
      int j = (i + 1) % fNv;
      good = AddFacet(new UQuadrangularFacet
                      (GetVertex(iz, j), GetVertex(iz, i),
                       GetVertex(iz + 1, i), GetVertex(iz + 1, j), UABSOLUTE));
      if (! good)
      {
        return false;
      }
    }
  }

  SetSolidClosed(true);

  return good;
}

//_____________________________________________________________________________

bool UExtrudedSolid::IsConvex() const
{
  // Get polygon convexity (polygon is convex if all vertex angles are < pi )

  for (int i = 0; i < fNv; ++i)
  {
    int j = (i + 1) % fNv;
    int k = (i + 2) % fNv;
    UVector2 v1 = fPolygon[i] - fPolygon[j];
    UVector2 v2 = fPolygon[k] - fPolygon[j];
    double dphi = v2.phi() - v1.phi();
    if (dphi < 0.)
    {
      dphi += 2.*UUtils::kPi;
    }

    if (dphi >= UUtils::kPi)
    {
      return false;
    }
  }

  return true;
}

//_____________________________________________________________________________

UGeometryType UExtrudedSolid::GetEntityType () const
{
  // Return entity type

  return fGeometryType;
}

//_____________________________________________________________________________

VUSolid* UExtrudedSolid::Clone() const
{
  return new UExtrudedSolid(*this);
}

//_____________________________________________________________________________

VUSolid::EnumInside UExtrudedSolid::Inside(const UVector3& p) const
{
  // Override the base class function  as it fails in case of concave polygon.
  // Project the point in the original polygon scale and check if it is inside
  // for each triangle.

  // Check first if outside extent
  //
  if (p.x() < GetMinXExtent() - VUSolid::fgTolerance * 0.5 ||
      p.x() > GetMaxXExtent() + VUSolid::fgTolerance * 0.5 ||
      p.y() < GetMinYExtent() - VUSolid::fgTolerance * 0.5 ||
      p.y() > GetMaxYExtent() + VUSolid::fgTolerance * 0.5 ||
      p.z() < GetMinZExtent() - VUSolid::fgTolerance * 0.5 ||
      p.z() > GetMaxZExtent() + VUSolid::fgTolerance * 0.5)
  {
    return eOutside;
  }

  // Project point p(z) to the polygon scale p0
  //
  UVector2 pscaled = ProjectPoint(p);

  // Check if on surface of polygon
  //
  for (int i = 0; i < fNv; ++i)
  {
    int j = (i + 1) % fNv;
    if (IsSameLineSegment(pscaled, fPolygon[i], fPolygon[j]))
    {
      return eSurface;
    }
  }

  // Now check if inside triangles
  //
  std::vector< std::vector<int> >::const_iterator it = fTriangles.begin();
  bool inside = false;
  do
  {
    if (IsPointInside(fPolygon[(*it)[0]], fPolygon[(*it)[1]],
                      fPolygon[(*it)[2]], pscaled))
    {
      inside = true;
    }
    ++it;
  }
  while ((inside == false) && (it != fTriangles.end()));

  if (inside)
  {
    // Check if on surface of z sides
    //
    if (std::fabs(p.z() - fZSections[0].fZ) < VUSolid::fgTolerance * 0.5 ||
        std::fabs(p.z() - fZSections[fNz - 1].fZ) < VUSolid::fgTolerance * 0.5)
    {
      return eSurface;
    }
    return eInside;
  }
  return eOutside;
}

//_____________________________________________________________________________
double UExtrudedSolid::DistanceToOut(const UVector3& p, const UVector3&  v, UVector3&       aNormalVector, bool&           aConvex, double) const
//double UExtrudedSolid::DistanceToOut (const UVector3 &p,
//                                       const UVector3 &v,
//                                       const bool calcNorm,
//                                             bool *validNorm,
//                                             UVector3 *n) const
{
  // Override the base class function to redefine validNorm
  // (the solid can be concave)

  double distOut =
    UTessellatedSolid::DistanceToOut(p, v, aNormalVector, aConvex);
  aConvex = fIsConvex;

  return distOut;
}


//_____________________________________________________________________________

double UExtrudedSolid::SafetyFromInside(const UVector3& p, bool) const
{
  // Override the overloaded base class function

  return UTessellatedSolid::SafetyFromInside(p);
}

//_____________________________________________________________________________

std::ostream& UExtrudedSolid::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid geometry type: " << fGeometryType  << std::endl;

  if (fIsConvex)
  {
    os << " Convex polygon; list of vertices:" << std::endl;
  }
  else
  {
    os << " Concave polygon; list of vertices:" << std::endl;
  }

  for (int i = 0; i < fNv; ++i)
  {
    os << std::setw(5) << "#" << i
       << "   vx = " << fPolygon[i].x << " mm"
       << "   vy = " << fPolygon[i].y << " mm" << std::endl;
  }

  os << " Sections:" << std::endl;
  for (int iz = 0; iz < fNz; ++iz)
  {
    os << "   z = "   << fZSections[iz].fZ          << " mm  "
       << "  x0= "    << fZSections[iz].fOffset.x << " mm  "
       << "  y0= "    << fZSections[iz].fOffset.y << " mm  "
       << "  scale= " << fZSections[iz].fScale << std::endl;
  }

  /*
    // Triangles (for debugging)
    os << std::endl;
    os << " Triangles:" << std::endl;
    os << " Triangle #   vertex1   vertex2   vertex3" << std::endl;

    int counter = 0;
    std::vector< std::vector<int> >::const_iterator it;
    for ( it = fTriangles.begin(); it != fTriangles.end(); it++ ) {
       std::vector<int> triangle = *it;
       os << std::setw(10) << counter++
          << std::setw(10) << triangle[0] << std::setw(10)  << triangle[1]  << std::setw(10)  << triangle[2]
          << std::endl;
    }
  */
  os.precision(oldprc);

  return os;
}
