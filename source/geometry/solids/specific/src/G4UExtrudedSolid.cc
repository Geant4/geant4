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
// Implementation of G4UExtrudedSolid wrapper class
//
// 17.11.17 G.Cosmo, CERN
// --------------------------------------------------------------------

#include "G4ExtrudedSolid.hh"
#include "G4UExtrudedSolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   const std::vector<G4TwoVector>& polygon,
                                   const std::vector<ZSection>&    zsections)
  : Base_t(name)  // General constructor
{
  unsigned int nVertices = polygon.size();
  unsigned int nSections = zsections.size();

  vecgeom::XtruVertex2* vertices = new vecgeom::XtruVertex2[nVertices];
  vecgeom::XtruSection* sections = new vecgeom::XtruSection[nSections];

  for (unsigned int i = 0; i < nVertices; ++i)
  {
    vertices[i].x = polygon[i].x();
    vertices[i].y = polygon[i].y();
  }
  for (unsigned int i = 0; i < nSections; ++i)
  {
    sections[i].fOrigin.Set(zsections[i].fOffset.x(),
                            zsections[i].fOffset.y(),
                            zsections[i].fZ);
    sections[i].fScale = zsections[i].fScale;
  }
  Base_t::Initialize(nVertices, vertices, nSections, sections);
  delete[] vertices;
  delete[] sections;
}


G4UExtrudedSolid::G4UExtrudedSolid(const G4String&          name,
                                   const std::vector<G4TwoVector>& polygon,
                                   G4double                 halfZ,
                                   const G4TwoVector& off1, G4double scale1,
                                   const G4TwoVector& off2, G4double scale2)
  : Base_t(name)  // Special constructor for 2 sections
{
  unsigned int nVertices = polygon.size();
  unsigned int nSections = 2;

  vecgeom::XtruVertex2* vertices = new vecgeom::XtruVertex2[nVertices];
  vecgeom::XtruSection* sections = new vecgeom::XtruSection[nSections];

  for (unsigned int i = 0; i < nVertices; ++i)
  {
    vertices[i].x = polygon[i].x();
    vertices[i].y = polygon[i].y();
  }
  sections[0].fOrigin.Set(off1.x(), off1.y(), -halfZ);
  sections[0].fScale = scale1;
  sections[1].fOrigin.Set(off2.x(), off2.y(), halfZ);
  sections[1].fScale = scale2;
  Base_t::Initialize(nVertices, vertices, nSections, sections);
  delete[] vertices;
  delete[] sections;
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UExtrudedSolid::G4UExtrudedSolid(__void__& a)
  : Base_t(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UExtrudedSolid::~G4UExtrudedSolid()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UExtrudedSolid::G4UExtrudedSolid(const G4UExtrudedSolid &source)
  : Base_t(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UExtrudedSolid&
G4UExtrudedSolid::operator=(const G4UExtrudedSolid &source)
{
  if (this == &source) return *this;
  
  Base_t::operator=( source );
  
  return *this;
}


//////////////////////////////////////////////////////////////////////////
//
// Accessors

G4int G4UExtrudedSolid::GetNofVertices() const
{
  return Base_t::GetNVertices();
}
G4TwoVector G4UExtrudedSolid::GetVertex(G4int i) const
{
  G4double xx, yy;
  Base_t::GetVertex(i, xx, yy);
  return G4TwoVector(xx, yy);
}
std::vector<G4TwoVector> G4UExtrudedSolid::GetPolygon() const
{
  std::vector<G4TwoVector> pol;
  for (unsigned int i = 0; i < Base_t::GetNVertices(); ++i)
  {
    pol.push_back(GetVertex(i));
  }
  return pol;
}
G4int G4UExtrudedSolid::GetNofZSections() const
{
  return Base_t::GetNSections();
}
G4UExtrudedSolid::ZSection G4UExtrudedSolid::GetZSection(G4int i) const
{
  vecgeom::XtruSection sect = Base_t::GetSection(i);
  return ZSection(sect.fOrigin[2],
                  G4TwoVector(sect.fOrigin[0], sect.fOrigin[1]),
                  sect.fScale);
}
std::vector<G4UExtrudedSolid::ZSection> G4UExtrudedSolid::GetZSections() const
{
  std::vector<G4UExtrudedSolid::ZSection> sections;
  for (unsigned int i = 0; i < Base_t::GetNSections(); ++i)
  {
    vecgeom::XtruSection sect = Base_t::GetSection(i);
    sections.push_back(ZSection(sect.fOrigin[2],
                                G4TwoVector(sect.fOrigin[0], sect.fOrigin[1]),
                                sect.fScale));
  }
  return sections;
}


///////////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UExtrudedSolid::BoundingLimits(G4ThreeVector& pMin,
                                      G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double xmin0 = kInfinity, xmax0 = -kInfinity;
  G4double ymin0 = kInfinity, ymax0 = -kInfinity;

  for (G4int i=0; i<GetNofVertices(); ++i)
  {
    G4TwoVector vertex = GetVertex(i);
    G4double x = vertex.x();
    if (x < xmin0) xmin0 = x;
    if (x > xmax0) xmax0 = x;
    G4double y = vertex.y();
    if (y < ymin0) ymin0 = y;
    if (y > ymax0) ymax0 = y;
  }

  G4double xmin = kInfinity, xmax = -kInfinity;
  G4double ymin = kInfinity, ymax = -kInfinity;

  G4int nsect = GetNofZSections();
  for (G4int i=0; i<nsect; ++i)
  {
    ZSection zsect = GetZSection(i);
    G4double dx    = zsect.fOffset.x();
    G4double dy    = zsect.fOffset.y();
    G4double scale = zsect.fScale;
    xmin = std::min(xmin,xmin0*scale+dx);
    xmax = std::max(xmax,xmax0*scale+dx);
    ymin = std::min(ymin,ymin0*scale+dy);
    ymax = std::max(ymax,ymax0*scale+dy);
  }

  G4double zmin = GetZSection(0).fZ;
  G4double zmax = GetZSection(nsect-1).fZ;

  pMin.set(xmin,ymin,zmin);
  pMax.set(xmax,ymax,zmax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UExtrudedSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    U3Vector vmin, vmax;
    Base_t::Extent(vmin,vmax);
    if (std::abs(pMin.x()-vmin.x()) > kCarTolerance ||
        std::abs(pMin.y()-vmin.y()) > kCarTolerance ||
        std::abs(pMin.z()-vmin.z()) > kCarTolerance ||
        std::abs(pMax.x()-vmax.x()) > kCarTolerance ||
        std::abs(pMax.y()-vmax.y()) > kCarTolerance ||
        std::abs(pMax.z()-vmax.z()) > kCarTolerance)
    {
      std::ostringstream message;
      message << "Inconsistency in bounding boxes for solid: "
              << GetName() << " !"
              << "\nBBox min: wrapper = " << pMin << " solid = " << vmin
              << "\nBBox max: wrapper = " << pMax << " solid = " << vmax;
      G4Exception("G4UExtrudedSolid::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}


//////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UExtrudedSolid::CalculateExtent(const EAxis pAxis,
                                  const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform,
                                        G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
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

  // To find the extent, the base polygon is subdivided in triangles.
  // The extent is calculated as cumulative extent of the parts
  // formed by extrusion of the triangles
  //
  G4TwoVectorList basePolygon = GetPolygon();
  G4TwoVectorList triangles;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  // triangulate the base polygon
  if (!G4GeomTools::TriangulatePolygon(basePolygon,triangles))
  {
    std::ostringstream message;
    message << "Triangulation of the base polygon has failed for solid: "
            << GetName() << " !"
            << "\nExtent has been calculated using boundary box";
    G4Exception("G4UExtrudedSolid::CalculateExtent()",
                "GeomMgt1002",JustWarning,message);
    return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }

  // allocate vector lists
  G4int nsect = GetNofZSections();
  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(nsect);
  for (G4int k=0; k<nsect; ++k) { polygons[k] = new G4ThreeVectorList(3); }

  // main loop along triangles
  pMin =  kInfinity;
  pMax = -kInfinity;
  G4int ntria = triangles.size()/3;
  for (G4int i=0; i<ntria; ++i)
  {
    G4int i3 = i*3;
    for (G4int k=0; k<nsect; ++k) // extrude triangle
    {
      ZSection zsect = GetZSection(k);
      G4double z     = zsect.fZ;
      G4double dx    = zsect.fOffset.x();
      G4double dy    = zsect.fOffset.y();
      G4double scale = zsect.fScale;

      G4ThreeVectorList* ptr = const_cast<G4ThreeVectorList*>(polygons[k]);
      G4ThreeVectorList::iterator iter = ptr->begin();
      G4double x0 = triangles[i3+0].x()*scale+dx;
      G4double y0 = triangles[i3+0].y()*scale+dy;
      iter->set(x0,y0,z);
      iter++;
      G4double x1 = triangles[i3+1].x()*scale+dx;
      G4double y1 = triangles[i3+1].y()*scale+dy;
      iter->set(x1,y1,z);
      iter++;
      G4double x2 = triangles[i3+2].x()*scale+dx;
      G4double y2 = triangles[i3+2].y()*scale+dy;
      iter->set(x2,y2,z);
    }

    // set sub-envelope and adjust extent
    G4double emin,emax;
    G4BoundingEnvelope benv(polygons);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) break; // max possible extent
  }
  // free memory
  for (G4int k=0; k<nsect; ++k) { delete polygons[k]; polygons[k]=0;}
  return (pMin < pMax);
}


///////////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron()
//
G4Polyhedron* G4UExtrudedSolid::CreatePolyhedron () const
{
  unsigned int nFacets = Base_t::GetStruct().fTslHelper.fFacets.size();
  unsigned int nVertices = Base_t::GetStruct().fTslHelper.fVertices.size();

  G4Polyhedron* polyhedron = new G4Polyhedron(nVertices, nFacets);

  // Copy vertices
  for (unsigned int i = 0; i < nVertices; ++i)
  {
    U3Vector v = Base_t::GetStruct().fTslHelper.fVertices[i];
    polyhedron->SetVertex(i+1, G4ThreeVector(v.x(), v.y(), v.z()));
  }

  // Copy facets
  for (unsigned int i = 0; i < nFacets; ++i)
  {
    // Facets are only triangular in VecGeom
    G4int i1 = Base_t::GetStruct().fTslHelper.fFacets[i]->fIndices[0] + 1;
    G4int i2 = Base_t::GetStruct().fTslHelper.fFacets[i]->fIndices[1] + 1;
    G4int i3 = Base_t::GetStruct().fTslHelper.fFacets[i]->fIndices[2] + 1;
    polyhedron->SetFacet(i+1, i1, i2, i3);
  }
  polyhedron->SetReferences();

  return polyhedron;
}

#endif  // G4GEOM_USE_USOLIDS
