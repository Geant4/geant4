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
// G4ExtrudedSolid implementation
//
// Author: Ivana Hrivnacova, IPN Orsay
//
// CHANGE HISTORY
// --------------
//
// 31.10.2017 E.Tcherniaev: added implementation for a non-convex
//            right prism
// 08.09.2017 E.Tcherniaev: added implementation for a convex
//            right prism
// 21.10.2016 E.Tcherniaev: reimplemented CalculateExtent(),
//            used G4GeomTools::PolygonArea() to calculate area,
//            replaced IsConvex() with G4GeomTools::IsConvex()
// 02.03.2016 E.Tcherniaev: added CheckPolygon() to remove
//            collinear and coincident points from polygon
// --------------------------------------------------------------------

#include "G4ExtrudedSolid.hh"

#if !defined(G4GEOM_USE_UEXTRUDEDSOLID)

#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( const G4String& pName,
                                  const std::vector<G4TwoVector>& polygon,
                                  const std::vector<ZSection>& zsections)
  : G4TessellatedSolid(pName),
    fNv(polygon.size()),
    fNz(zsections.size()),
    fIsConvex(false),
    fGeometryType("G4ExtrudedSolid"),
    fSolidType(0)
{
  // General constructor

  // First check input parameters

  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of vertices in polygon < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  if (fNz < 2)
  {
    std::ostringstream message;
    message << "Number of z-sides < 2 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  for ( std::size_t i=0; i<fNz-1; ++i ) 
  {
    if ( zsections[i].fZ > zsections[i+1].fZ ) 
    {
      std::ostringstream message;
      message << "Z-sections have to be ordered by z value (z0 < z1 < z2...) - "
              << pName;
      G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                  FatalErrorInArgument, message);
    }
    if ( std::fabs( zsections[i+1].fZ - zsections[i].fZ ) < kCarToleranceHalf )
    {
      std::ostringstream message;
      message << "Z-sections with the same z position are not supported - "
              << pName;
      G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0001",
                  FatalException, message);
    }
  }  
  
  // Copy polygon
  //
  fPolygon = polygon;

  // Remove collinear and coincident vertices, if any
  //
  std::vector<G4int> removedVertices;
  G4GeomTools::RemoveRedundantVertices(fPolygon,removedVertices,
                                       2*kCarTolerance);
  if (removedVertices.size() != 0)
  {
    std::size_t nremoved = removedVertices.size();
    std::ostringstream message;
    message << "The following "<< nremoved 
            << " vertices have been removed from polygon in " << pName 
            << "\nas collinear or coincident with other vertices: "
            << removedVertices[0];
    for (std::size_t i=1; i<nremoved; ++i) message << ", " << removedVertices[i];
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
                JustWarning, message);
  }

  fNv = fPolygon.size();
  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of vertices in polygon after removal < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)
  //
  if (G4GeomTools::PolygonArea(fPolygon) > 0.)
  {   
    // Polygon vertices are defined anti-clockwise, we revert them
    // G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
    //            JustWarning, 
    //            "Polygon vertices defined anti-clockwise, reverting polygon");
    std::reverse(fPolygon.begin(),fPolygon.end());
  }

  // Copy z-sections
  //
  fZSections = zsections;

  G4bool result = MakeFacets();
  if (!result)
  {   
    std::ostringstream message;
    message << "Making facets failed - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0003",
                FatalException, message);
  }
  fIsConvex = G4GeomTools::IsConvex(fPolygon);
  
  ComputeProjectionParameters();

  // Check if the solid is a right prism, if so then set lateral planes
  //
  if ((fNz == 2)
      && (fZSections[0].fScale == 1) && (fZSections[1].fScale == 1)
      && (fZSections[0].fOffset == G4TwoVector(0,0))
      && (fZSections[1].fOffset == G4TwoVector(0,0)))
  {
    fSolidType = (fIsConvex) ? 1 : 2; // 1 - convex, 2 - non-convex right prism
    ComputeLateralPlanes();
  }
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( const G4String& pName,
                                  const std::vector<G4TwoVector>& polygon,
                                        G4double dz,
                                  const G4TwoVector& off1, G4double scale1,
                                  const G4TwoVector& off2, G4double scale2 )
  : G4TessellatedSolid(pName),
    fNv(polygon.size()),
    fNz(2),
    fGeometryType("G4ExtrudedSolid")
{
  // Special constructor for solid with 2 z-sections

  // First check input parameters
  //
  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of vertices in polygon < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  // Copy polygon
  //
  fPolygon = polygon;

  // Remove collinear and coincident vertices, if any
  //
  std::vector<G4int> removedVertices;
  G4GeomTools::RemoveRedundantVertices(fPolygon,removedVertices,
                                       2*kCarTolerance);
  if (removedVertices.size() != 0)
  {
    std::size_t nremoved = removedVertices.size();
    std::ostringstream message;
    message << "The following "<< nremoved 
            << " vertices have been removed from polygon in " << pName 
            << "\nas collinear or coincident with other vertices: "
            << removedVertices[0];
    for (std::size_t i=1; i<nremoved; ++i) message << ", " << removedVertices[i];
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
                JustWarning, message);
  }

  fNv = fPolygon.size();
  if (fNv < 3)
  {
    std::ostringstream message;
    message << "Number of vertices in polygon after removal < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)
  //  
  if (G4GeomTools::PolygonArea(fPolygon) > 0.)
  {
    // Polygon vertices are defined anti-clockwise, we revert them
    // G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
    //            JustWarning, 
    //            "Polygon vertices defined anti-clockwise, reverting polygon");
    std::reverse(fPolygon.begin(),fPolygon.end());
  }
  
  // Copy z-sections
  //
  fZSections.push_back(ZSection(-dz, off1, scale1));
  fZSections.push_back(ZSection( dz, off2, scale2));
    
  G4bool result = MakeFacets();
  if (!result)
  {   
    std::ostringstream message;
    message << "Making facets failed - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0003",
                FatalException, message);
  }
  fIsConvex = G4GeomTools::IsConvex(fPolygon);

  ComputeProjectionParameters();

  // Check if the solid is a right prism, if so then set lateral planes
  //
  if ((scale1 == 1) && (scale2 == 1)
      && (off1 == G4TwoVector(0,0)) && (off2 == G4TwoVector(0,0)))
  {
    fSolidType = (fIsConvex) ? 1 : 2; // 1 - convex, 2 - non-convex right prism
    ComputeLateralPlanes();
  }
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( __void__& a )
  : G4TessellatedSolid(a), fNv(0), fNz(0),
    fGeometryType("G4ExtrudedSolid")
{
  // Fake default constructor - sets only member data and allocates memory
  //                            for usage restricted to object persistency.
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid(const G4ExtrudedSolid& rhs)
  : G4TessellatedSolid(rhs), fNv(rhs.fNv), fNz(rhs.fNz),
    fPolygon(rhs.fPolygon), fZSections(rhs.fZSections),
    fTriangles(rhs.fTriangles), fIsConvex(rhs.fIsConvex),
    fGeometryType(rhs.fGeometryType),
    fSolidType(rhs.fSolidType), fPlanes(rhs.fPlanes),
    fLines(rhs.fLines), fLengths(rhs.fLengths),
    fKScales(rhs.fKScales), fScale0s(rhs.fScale0s),
    fKOffsets(rhs.fKOffsets), fOffset0s(rhs.fOffset0s)
{
}

//_____________________________________________________________________________

G4ExtrudedSolid& G4ExtrudedSolid::operator = (const G4ExtrudedSolid& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4TessellatedSolid::operator=(rhs);

   // Copy data
   //
   fNv = rhs.fNv; fNz = rhs.fNz;
   fPolygon = rhs.fPolygon; fZSections = rhs.fZSections;
   fTriangles = rhs.fTriangles; fIsConvex = rhs.fIsConvex;
   fGeometryType = rhs.fGeometryType;
   fSolidType = rhs.fSolidType; fPlanes = rhs.fPlanes;
   fLines = rhs.fLines; fLengths = rhs.fLengths;
   fKScales = rhs.fKScales; fScale0s = rhs.fScale0s;
   fKOffsets = rhs.fKOffsets; fOffset0s = rhs.fOffset0s;

   return *this;
}

//_____________________________________________________________________________

G4ExtrudedSolid::~G4ExtrudedSolid()
{
  // Destructor
}

//_____________________________________________________________________________

void G4ExtrudedSolid::ComputeProjectionParameters()
{
  // Compute parameters for point projections p(z)
  // to the polygon scale & offset:
  // scale(z) = k*z + scale0
  // offset(z) = l*z + offset0
  // p(z) = scale(z)*p0 + offset(z)
  // p0 = (p(z) - offset(z))/scale(z);
  //

  for (std::size_t iz=0; iz<fNz-1; ++iz)
  {
    G4double z1      = fZSections[iz].fZ;
    G4double z2      = fZSections[iz+1].fZ;
    G4double scale1  = fZSections[iz].fScale;
    G4double scale2  = fZSections[iz+1].fScale;
    G4TwoVector off1 = fZSections[iz].fOffset;
    G4TwoVector off2 = fZSections[iz+1].fOffset;

    G4double kscale = (scale2 - scale1)/(z2 - z1);
    G4double scale0 =  scale2 - kscale*(z2 - z1)/2.0;
    G4TwoVector koff = (off2 - off1)/(z2 - z1);
    G4TwoVector off0 =  off2 - koff*(z2 - z1)/2.0;

    fKScales.push_back(kscale);
    fScale0s.push_back(scale0);
    fKOffsets.push_back(koff);
    fOffset0s.push_back(off0);
  }
}

//_____________________________________________________________________________

void G4ExtrudedSolid::ComputeLateralPlanes()
{
  // Compute lateral planes: a*x + b*y + c*z + d = 0
  //
  std::size_t Nv = fPolygon.size();
  fPlanes.resize(Nv);
  for (std::size_t i=0, k=Nv-1; i<Nv; k=i++)
  {
    G4TwoVector norm = (fPolygon[i] - fPolygon[k]).unit();
    fPlanes[i].a = -norm.y();
    fPlanes[i].b =  norm.x();
    fPlanes[i].c =  0;
    fPlanes[i].d =  norm.y()*fPolygon[i].x() - norm.x()*fPolygon[i].y();
  }

  // Compute edge equations: x = k*y + m
  // and edge lengths
  //
  fLines.resize(Nv);
  fLengths.resize(Nv);
  for (std::size_t i=0, k=Nv-1; i<Nv; k=i++)
  {
    if (fPolygon[k].y() == fPolygon[i].y())
    {
      fLines[i].k = 0;
      fLines[i].m = fPolygon[i].x();
    }
    else
    {
      G4double ctg = (fPolygon[k].x()-fPolygon[i].x())/(fPolygon[k].y()-fPolygon[i].y());
      fLines[i].k = ctg;
      fLines[i].m = fPolygon[i].x() - ctg*fPolygon[i].y();
    }
    fLengths[i]  =  (fPolygon[i] - fPolygon[k]).mag();
  }
}

//_____________________________________________________________________________

G4ThreeVector G4ExtrudedSolid::GetVertex(G4int iz, G4int ind) const
{
  // Shift and scale vertices

  return G4ThreeVector( fPolygon[ind].x() * fZSections[iz].fScale
                      + fZSections[iz].fOffset.x(), 
                        fPolygon[ind].y() * fZSections[iz].fScale
                      + fZSections[iz].fOffset.y(), fZSections[iz].fZ);
}       

//_____________________________________________________________________________

G4TwoVector G4ExtrudedSolid::ProjectPoint(const G4ThreeVector& point) const
{
  // Project point in the polygon scale
  // scale(z) = k*z + scale0
  // offset(z) = l*z + offset0
  // p(z) = scale(z)*p0 + offset(z)  
  // p0 = (p(z) - offset(z))/scale(z);
  
  // Select projection (z-segment of the solid) according to p.z()
  //
  std::size_t iz = 0;
  while ( point.z() > fZSections[iz+1].fZ && iz < fNz-2 ) { ++iz; }
      // Loop checking, 13.08.2015, G.Cosmo
  
  G4double z0 = ( fZSections[iz+1].fZ + fZSections[iz].fZ )/2.0;
  G4TwoVector p2(point.x(), point.y());
  G4double pscale  = fKScales[iz]*(point.z()-z0) + fScale0s[iz];
  G4TwoVector poffset = fKOffsets[iz]*(point.z()-z0) + fOffset0s[iz];
  
  // G4cout << point << " projected to " 
  //        << iz << "-th z-segment polygon as "
  //        << (p2 - poffset)/pscale << G4endl;

  // pscale is always >0 as it is an interpolation between two
  // positive scale values
  //
  return (p2 - poffset)/pscale;
}  

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsSameLine(const G4TwoVector& p,
                                   const G4TwoVector& l1,
                                   const G4TwoVector& l2) const
{
  // Return true if p is on the line through l1, l2 

  if ( l1.x() == l2.x() )
  {
    return std::fabs(p.x() - l1.x()) < kCarToleranceHalf; 
  }
   G4double slope= ((l2.y() - l1.y())/(l2.x() - l1.x())); 
   G4double predy= l1.y() +  slope *(p.x() - l1.x()); 
   G4double dy= p.y() - predy; 

   // Calculate perpendicular distance
   //
   // G4double perpD= std::fabs(dy) / std::sqrt( 1 + slope * slope ); 
   // G4bool   simpleComp= (perpD<kCarToleranceHalf); 

   // Check perpendicular distance vs tolerance 'directly'
   //
   G4bool squareComp = (dy*dy < (1+slope*slope)
                     * kCarToleranceHalf * kCarToleranceHalf); 

   // return  simpleComp; 
   return squareComp;
}                    

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsSameLineSegment(const G4TwoVector& p,  
                                          const G4TwoVector& l1,
                                          const G4TwoVector& l2) const
{
  // Return true if p is on the line through l1, l2 and lies between
  // l1 and l2 

  if ( p.x() < std::min(l1.x(), l2.x()) - kCarToleranceHalf || 
       p.x() > std::max(l1.x(), l2.x()) + kCarToleranceHalf ||
       p.y() < std::min(l1.y(), l2.y()) - kCarToleranceHalf || 
       p.y() > std::max(l1.y(), l2.y()) + kCarToleranceHalf )
  {
    return false;
  }

  return IsSameLine(p, l1, l2);
}

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsSameSide(const G4TwoVector& p1,
                                   const G4TwoVector& p2,
                                   const G4TwoVector& l1,
                                   const G4TwoVector& l2) const
{
  // Return true if p1 and p2 are on the same side of the line through l1, l2 

  return   ( (p1.x() - l1.x()) * (l2.y() - l1.y())
         - (l2.x() - l1.x()) * (p1.y() - l1.y()) )
         * ( (p2.x() - l1.x()) * (l2.y() - l1.y())
         - (l2.x() - l1.x()) * (p2.y() - l1.y()) ) > 0;
}       

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsPointInside(const G4TwoVector& a,
                                      const G4TwoVector& b,
                                      const G4TwoVector& c,
                                      const G4TwoVector& p) const
{
  // Return true if p is inside of triangle abc or on its edges, 
  // else returns false 

  // Check extent first
  //
  if ( ( p.x() < a.x() && p.x() < b.x() && p.x() < c.x() ) || 
       ( p.x() > a.x() && p.x() > b.x() && p.x() > c.x() ) || 
       ( p.y() < a.y() && p.y() < b.y() && p.y() < c.y() ) || 
       ( p.y() > a.y() && p.y() > b.y() && p.y() > c.y() ) ) return false;
  
  G4bool inside 
    = IsSameSide(p, a, b, c)
      && IsSameSide(p, b, a, c)
      && IsSameSide(p, c, a, b);

  G4bool onEdge
    = IsSameLineSegment(p, a, b)
      || IsSameLineSegment(p, b, c)
      || IsSameLineSegment(p, c, a);
      
  return inside || onEdge;    
}     

//_____________________________________________________________________________

G4double 
G4ExtrudedSolid::GetAngle(const G4TwoVector& po,
                          const G4TwoVector& pa,
                          const G4TwoVector& pb) const
{
  // Return the angle of the vertex in po

  G4TwoVector t1 = pa - po;
  G4TwoVector t2 = pb - po;
  
  G4double result = (std::atan2(t1.y(), t1.x()) - std::atan2(t2.y(), t2.x()));

  if ( result < 0 ) result += 2*pi;

  return result;
}

//_____________________________________________________________________________

G4VFacet*
G4ExtrudedSolid::MakeDownFacet(G4int ind1, G4int ind2, G4int ind3) const
{
  // Create a triangular facet from the polygon points given by indices
  // forming the down side ( the normal goes in -z)

  std::vector<G4ThreeVector> vertices;
  vertices.push_back(GetVertex(0, ind1));
  vertices.push_back(GetVertex(0, ind2));
  vertices.push_back(GetVertex(0, ind3));
  
  // first vertex most left
  //
  G4ThreeVector cross 
    = (vertices[1]-vertices[0]).cross(vertices[2]-vertices[1]);
  
  if ( cross.z() > 0.0 )
  {
    // vertices ordered clock wise has to be reordered

    // G4cout << "G4ExtrudedSolid::MakeDownFacet: reordering vertices " 
    //        << ind1 << ", " << ind2 << ", " << ind3 << G4endl; 

    G4ThreeVector tmp = vertices[1];
    vertices[1] = vertices[2];
    vertices[2] = tmp;
  }
  
  return new G4TriangularFacet(vertices[0], vertices[1],
                               vertices[2], ABSOLUTE);
}      

//_____________________________________________________________________________

G4VFacet*
G4ExtrudedSolid::MakeUpFacet(G4int ind1, G4int ind2, G4int ind3) const      
{
  // Creates a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  std::vector<G4ThreeVector> vertices;
  vertices.push_back(GetVertex((G4int)fNz-1, ind1));
  vertices.push_back(GetVertex((G4int)fNz-1, ind2));
  vertices.push_back(GetVertex((G4int)fNz-1, ind3));
  
  // first vertex most left
  //
  G4ThreeVector cross 
    = (vertices[1]-vertices[0]).cross(vertices[2]-vertices[1]);
  
  if ( cross.z() < 0.0 )
  {
    // vertices ordered clock wise has to be reordered

    // G4cout << "G4ExtrudedSolid::MakeUpFacet: reordering vertices " 
    //        << ind1 << ", " << ind2 << ", " << ind3 << G4endl; 

    G4ThreeVector tmp = vertices[1];
    vertices[1] = vertices[2];
    vertices[2] = tmp;
  }
  
  return new G4TriangularFacet(vertices[0], vertices[1],
                               vertices[2], ABSOLUTE);
}      

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::AddGeneralPolygonFacets()
{
  // Decompose polygonal sides in triangular facets

  typedef std::pair < G4TwoVector, G4int > Vertex;

  static const G4double kAngTolerance =
    G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  // Fill one more vector
  //
  std::vector< Vertex > verticesToBeDone;
  for ( G4int i=0; i<(G4int)fNv; ++i )
  {
    verticesToBeDone.push_back(Vertex(fPolygon[i], i));
  }
  std::vector< Vertex > ears;
  
  std::vector< Vertex >::iterator c1 = verticesToBeDone.begin();
  std::vector< Vertex >::iterator c2 = c1+1;  
  std::vector< Vertex >::iterator c3 = c1+2;  
  while ( verticesToBeDone.size()>2 )    // Loop checking, 13.08.2015, G.Cosmo
  {

    // G4cout << "Looking at triangle : "
    //         << c1->second << "  " << c2->second
    //        << "  " << c3->second << G4endl;  
    //G4cout << "Looking at triangle : "
    //        << c1->first << "  " << c2->first
    //        << "  " << c3->first << G4endl;  

    // skip concave vertices
    //
    G4double angle = GetAngle(c2->first, c3->first, c1->first);
   
    //G4cout << "angle " << angle  << G4endl;

    std::size_t counter = 0;
    while ( angle >= (pi-kAngTolerance) )  // Loop checking, 13.08.2015, G.Cosmo
    {
      // G4cout << "Skipping concave vertex " << c2->second << G4endl;

      // try next three consecutive vertices
      //
      c1 = c2;
      c2 = c3;
      ++c3; 
      if ( c3 == verticesToBeDone.end() ) { c3 = verticesToBeDone.begin(); }

      //G4cout << "Looking at triangle : "
      //      << c1->first << "  " << c2->first
      //        << "  " << c3->first << G4endl; 
      
      angle = GetAngle(c2->first, c3->first, c1->first); 
      //G4cout << "angle " << angle  << G4endl;
      
      ++counter;
      
      if ( counter > fNv )
      {
        G4Exception("G4ExtrudedSolid::AddGeneralPolygonFacets",
                    "GeomSolids0003", FatalException,
                    "Triangularisation has failed.");
        break;
      }  
    }

    G4bool good = true;
    for ( auto it=verticesToBeDone.cbegin(); it!=verticesToBeDone.cend(); ++it )
    {
      // skip vertices of tested triangle
      //
      if ( it == c1 || it == c2 || it == c3 ) { continue; }

      if ( IsPointInside(c1->first, c2->first, c3->first, it->first) )
      {
        // G4cout << "Point " << it->second << " is inside" << G4endl;
        good = false;

        // try next three consecutive vertices
        //
        c1 = c2;
        c2 = c3;
        ++c3; 
        if ( c3 == verticesToBeDone.end() ) { c3 = verticesToBeDone.begin(); }
        break;
      }
      // else 
      //   { G4cout << "Point " << it->second << " is outside" << G4endl; }
    }
    if ( good )
    {
      // all points are outside triangle, we can make a facet

      // G4cout << "Found triangle : "
      //        << c1->second << "  " << c2->second
      //        << "  " << c3->second << G4endl;  

      G4bool result;
      result = AddFacet( MakeDownFacet(c1->second, c2->second, c3->second) );
      if ( ! result ) { return false; }

      result = AddFacet( MakeUpFacet(c1->second, c2->second, c3->second) );
      if ( ! result ) { return false; }

      std::vector<G4int> triangle(3);
      triangle[0] = c1->second;
      triangle[1] = c2->second;
      triangle[2] = c3->second;
      fTriangles.push_back(triangle);

      // remove the ear point from verticesToBeDone
      //
      verticesToBeDone.erase(c2);
      c1 = verticesToBeDone.begin();
      c2 = c1+1;  
      c3 = c1+2;  
    } 
  }
  return true;
}

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::MakeFacets()
{
  // Define facets

  G4bool good;
  
  // Decomposition of polygonal sides in the facets
  //
  if ( fNv == 3 )
  {
    good = AddFacet( new G4TriangularFacet( GetVertex(0, 0), GetVertex(0, 1),
                                            GetVertex(0, 2), ABSOLUTE) );
    if ( ! good ) { return false; }

    good = AddFacet( new G4TriangularFacet( GetVertex((G4int)fNz-1, 2),
                                            GetVertex((G4int)fNz-1, 1),
                                            GetVertex((G4int)fNz-1, 0),
                                            ABSOLUTE) );
    if ( ! good ) { return false; }
    
    std::vector<G4int> triangle(3);
    triangle[0] = 0;
    triangle[1] = 1;
    triangle[2] = 2;
    fTriangles.push_back(triangle);
  }
  
  else if ( fNv == 4 )
  {
    good = AddFacet( new G4QuadrangularFacet( GetVertex(0, 0),GetVertex(0, 1),
                                              GetVertex(0, 2),GetVertex(0, 3),
                                              ABSOLUTE) );
    if ( ! good ) { return false; }

    good = AddFacet( new G4QuadrangularFacet( GetVertex((G4int)fNz-1, 3),
                                              GetVertex((G4int)fNz-1, 2), 
                                              GetVertex((G4int)fNz-1, 1),
                                              GetVertex((G4int)fNz-1, 0),
                                              ABSOLUTE) );
    if ( ! good ) { return false; }

    std::vector<G4int> triangle1(3);
    triangle1[0] = 0;
    triangle1[1] = 1;
    triangle1[2] = 2;
    fTriangles.push_back(triangle1);

    std::vector<G4int> triangle2(3);
    triangle2[0] = 0;
    triangle2[1] = 2;
    triangle2[2] = 3;
    fTriangles.push_back(triangle2);
  }  
  else
  {
    good = AddGeneralPolygonFacets();
    if ( ! good ) { return false; }
  }
    
  // The quadrangular sides
  //
  for ( G4int iz = 0; iz < (G4int)fNz-1; ++iz ) 
  {
    for ( G4int i = 0; i < (G4int)fNv; ++i )
    {
      G4int j = (i+1) % fNv;
      good = AddFacet( new G4QuadrangularFacet
                        ( GetVertex(iz, j), GetVertex(iz, i), 
                          GetVertex(iz+1, i), GetVertex(iz+1, j), ABSOLUTE) );
      if ( ! good ) { return false; }
    }
  }  

  SetSolidClosed(true);

  return good;
}

//_____________________________________________________________________________

G4GeometryType G4ExtrudedSolid::GetEntityType () const
{
  // Return entity type

  return fGeometryType;
}

//_____________________________________________________________________________

G4VSolid* G4ExtrudedSolid::Clone() const
{
  return new G4ExtrudedSolid(*this);
}

//_____________________________________________________________________________

EInside G4ExtrudedSolid::Inside(const G4ThreeVector &p) const
{
  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      G4double dist = std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
      if (dist > kCarToleranceHalf)  { return kOutside; }

      std::size_t np = fPlanes.size();
      for (std::size_t i=0; i<np; ++i)
      {
        G4double dd = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
        if (dd > dist)  { dist = dd; }
      }
      if (dist > kCarToleranceHalf)  { return kOutside; }
      return (dist > -kCarToleranceHalf) ? kSurface : kInside;
    }
    case 2: // non-convex right prism
    {
      G4double distz = std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
      if (distz > kCarToleranceHalf)  { return kOutside; }

      G4bool in = PointInPolygon(p);
      if (distz > -kCarToleranceHalf && in) { return kSurface; }

      G4double dd = DistanceToPolygonSqr(p) - kCarToleranceHalf*kCarToleranceHalf;
      if (in)
      {
        return (dd >= 0) ? kInside : kSurface;
      }
      else
      {
        return (dd > 0) ? kOutside : kSurface;
      }
    }
  }

  // Override the base class function  as it fails in case of concave polygon.
  // Project the point in the original polygon scale and check if it is inside
  // for each triangle.

  // Check first if outside extent
  //
  if ( p.x() < GetMinXExtent() - kCarToleranceHalf ||
       p.x() > GetMaxXExtent() + kCarToleranceHalf ||
       p.y() < GetMinYExtent() - kCarToleranceHalf ||
       p.y() > GetMaxYExtent() + kCarToleranceHalf ||
       p.z() < GetMinZExtent() - kCarToleranceHalf ||
       p.z() > GetMaxZExtent() + kCarToleranceHalf )
  {
    // G4cout << "G4ExtrudedSolid::Outside extent: " << p << G4endl;
    return kOutside;
  }

  // Project point p(z) to the polygon scale p0
  //
  G4TwoVector pscaled = ProjectPoint(p);
  
  // Check if on surface of polygon
  //
  for ( G4int i=0; i<(G4int)fNv; ++i )
  {
    G4int j = (i+1) % fNv;
    if ( IsSameLineSegment(pscaled, fPolygon[i], fPolygon[j]) )
    {
      // G4cout << "G4ExtrudedSolid::Inside return Surface (on polygon) "
      //        << G4endl;

      return kSurface;
    }
  }

  // Now check if inside triangles
  //
  auto it = fTriangles.cbegin();
  G4bool inside = false;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    if ( IsPointInside(fPolygon[(*it)[0]], fPolygon[(*it)[1]],
                       fPolygon[(*it)[2]], pscaled) )  { inside = true; }
    ++it;
  } while ( (inside == false) && (it != fTriangles.cend()) );

  if ( inside )
  {
    // Check if on surface of z sides
    //
    if ( std::fabs( p.z() - fZSections[0].fZ ) < kCarToleranceHalf ||
         std::fabs( p.z() - fZSections[fNz-1].fZ ) < kCarToleranceHalf )
    {
      // G4cout << "G4ExtrudedSolid::Inside return Surface (on z side)"
      //        << G4endl;

      return kSurface;
    }

    // G4cout << "G4ExtrudedSolid::Inside return Inside" << G4endl;

    return kInside;
  }

  // G4cout << "G4ExtrudedSolid::Inside return Outside " << G4endl;

  return kOutside;
}

//_____________________________________________________________________________

G4ThreeVector G4ExtrudedSolid::SurfaceNormal(const G4ThreeVector& p) const
{
  G4int nsurf = 0;
  G4double nx = 0., ny = 0., nz = 0.;
  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      if (std::abs(p.z() - fZSections[0].fZ) <= kCarToleranceHalf)
      {
        nz = -1; ++nsurf;
      }
      if (std::abs(p.z() - fZSections[1].fZ) <= kCarToleranceHalf)
      {
        nz =  1; ++nsurf;
      }
      for (std::size_t i=0; i<fNv; ++i)
      {
        G4double dd = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
        if (std::abs(dd) > kCarToleranceHalf) continue;
        nx += fPlanes[i].a;
        ny += fPlanes[i].b;
        ++nsurf;
      }
      break;
    }
    case 2: // non-convex right prism
    {
      if (std::abs(p.z() - fZSections[0].fZ) <= kCarToleranceHalf)
      {
        nz = -1; ++nsurf;
      }
      if (std::abs(p.z() - fZSections[1].fZ) <= kCarToleranceHalf)
      {
        nz =  1; ++nsurf;
      }

      G4double sqrCarToleranceHalf = kCarToleranceHalf*kCarToleranceHalf;
      for (std::size_t i=0, k=fNv-1; i<fNv; k=i++)
      {
        G4double ix = p.x() - fPolygon[i].x();
        G4double iy = p.y() - fPolygon[i].y();
        G4double u  = fPlanes[i].a*iy - fPlanes[i].b*ix;
        if (u < 0)
        {
          if (ix*ix + iy*iy > sqrCarToleranceHalf) continue;
        }
        else if (u > fLengths[i])
        {
          G4double kx = p.x() - fPolygon[k].x();
          G4double ky = p.y() - fPolygon[k].y();
          if (kx*kx + ky*ky > sqrCarToleranceHalf) continue;
        }
        else
        {
          G4double dd = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
          if (dd*dd > sqrCarToleranceHalf) continue;
        }
        nx += fPlanes[i].a;
        ny += fPlanes[i].b;
        ++nsurf;
      }
      break;
    }
    default:
    {
      return G4TessellatedSolid::SurfaceNormal(p);
    }
  }

  // Return normal (right prism)
  //
  if (nsurf == 1)
  {
    return G4ThreeVector(nx,ny,nz);
  }
  else if (nsurf != 0) // edge or corner
  {
    return G4ThreeVector(nx,ny,nz).unit();
  }
  else
  {
    // Point is not on the surface, compute approximate normal
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
    G4cout.precision(oldprc) ;
    G4Exception("G4TesselatedSolid::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

//_____________________________________________________________________________

G4ThreeVector G4ExtrudedSolid::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  // This method is valid only for right prisms and
  // normally should not be called

  if (fSolidType == 1 || fSolidType == 2)
  {
    // Find distances to z-planes
    //
    G4double dz0 = fZSections[0].fZ - p.z();
    G4double dz1 = p.z() - fZSections[1].fZ;
    G4double ddz0 = dz0*dz0;
    G4double ddz1 = dz1*dz1;

    // Find nearest lateral side and distance to it
    //
    std::size_t iside = 0;
    G4double dd = DBL_MAX;
    for (std::size_t i=0, k=fNv-1; i<fNv; k=i++)
    {
      G4double ix = p.x() - fPolygon[i].x();
      G4double iy = p.y() - fPolygon[i].y();
      G4double u  = fPlanes[i].a*iy - fPlanes[i].b*ix;
      if (u < 0)
      {
        G4double tmp = ix*ix + iy*iy;
        if (tmp < dd) { dd = tmp; iside = i; }
      }
      else if (u > fLengths[i])
      {
        G4double kx = p.x() - fPolygon[k].x();
        G4double ky = p.y() - fPolygon[k].y();
        G4double tmp = kx*kx + ky*ky;
        if (tmp < dd)  { dd = tmp; iside = i; }
      }
      else
      {
        G4double tmp = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
        tmp *= tmp;
        if (tmp < dd)  { dd = tmp; iside = i; }
      }
    }

    // Find region
    //
    //  3  |   1   |  3
    // ----+-------+----
    //  2  |   0   |  2
    // ----+-------+----
    //  3  |   1   |  3
    //
    G4int iregion = 0;
    if (std::max(dz0,dz1) > 0) iregion = 1;

    G4bool in = PointInPolygon(p);
    if (!in) iregion += 2;

    // Return normal
    //
    switch (iregion)
    {
      case 0:
      {
        if (ddz0 <= ddz1 && ddz0 <= dd) return G4ThreeVector(0, 0,-1);
        if (ddz1 <= ddz0 && ddz1 <= dd) return G4ThreeVector(0, 0, 1);
        return G4ThreeVector(fPlanes[iside].a, fPlanes[iside].b, 0);
      }
      case 1:
      {
        return G4ThreeVector(0, 0, (dz0 > dz1) ? -1 : 1);
      }
      case 2:
      {
        return G4ThreeVector(fPlanes[iside].a, fPlanes[iside].b, 0);
      }
      case 3:
      {
        G4double dzmax = std::max(dz0,dz1);
        if (dzmax*dzmax > dd) return G4ThreeVector(0,0,(dz0 > dz1) ? -1 : 1);
        return G4ThreeVector(fPlanes[iside].a,fPlanes[iside].b, 0);
      }
    }
  }
  return G4ThreeVector(0,0,0);
}

//_____________________________________________________________________________

G4double G4ExtrudedSolid::DistanceToIn(const G4ThreeVector& p,
                                       const G4ThreeVector& v) const
{
  G4double z0 = fZSections[0].fZ;
  G4double z1 = fZSections[fNz-1].fZ;
  if ((p.z() <= z0 + kCarToleranceHalf) && v.z() <= 0) return kInfinity;
  if ((p.z() >= z1 - kCarToleranceHalf) && v.z() >= 0) return kInfinity;

  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      // Intersection with Z planes
      //
      G4double dz = (z1 - z0)*0.5;
      G4double pz = p.z() - dz - z0;

      G4double invz = (v.z() == 0) ? DBL_MAX : -1./v.z();
      G4double ddz = (invz < 0) ? dz : -dz;
      G4double tzmin = (pz + ddz)*invz;
      G4double tzmax = (pz - ddz)*invz;

      // Intersection with lateral planes
      //
      std::size_t np = fPlanes.size();
      G4double txmin = tzmin, txmax = tzmax;
      for (std::size_t i=0; i<np; ++i)
      { 
        G4double cosa = fPlanes[i].a*v.x()+fPlanes[i].b*v.y();
        G4double dist = fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].d;
        if (dist >= -kCarToleranceHalf)
        {
          if (cosa >= 0)  { return kInfinity; }
          G4double tmp  = -dist/cosa;
          if (txmin < tmp)  { txmin = tmp; }
        }
        else if (cosa > 0)
        {
          G4double tmp  = -dist/cosa;
          if (txmax > tmp)  { txmax = tmp; }
        } 
      }

      // Find distance
      //
      G4double tmin = txmin, tmax = txmax;
      if (tmax <= tmin + kCarToleranceHalf)   // touch or no hit
      {
        return kInfinity;
      }
      return (tmin < kCarToleranceHalf) ? 0. : tmin;
    }
    case 2: // non-convex right prism
    {
    }
  }
  return G4TessellatedSolid::DistanceToIn(p,v);
}

//_____________________________________________________________________________

G4double G4ExtrudedSolid::DistanceToIn (const G4ThreeVector& p) const
{
  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      G4double dist = std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
      std::size_t np = fPlanes.size();
      for (std::size_t i=0; i<np; ++i)
      {
        G4double dd = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
        if (dd > dist) dist = dd;
      }
      return (dist > 0) ? dist : 0.;
    }
    case 2: // non-convex right prism
    {
      G4bool in = PointInPolygon(p);
      if (in)
      {
        G4double distz= std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
        return (distz > 0) ? distz : 0;
      }
      else
      {
        G4double distz= std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
        G4double dd = DistanceToPolygonSqr(p);
        if (distz > 0) dd += distz*distz;
        return std::sqrt(dd);
      }
    }
  }

  // General case: use tessellated solid
  return G4TessellatedSolid::DistanceToIn(p);
}

//_____________________________________________________________________________

G4double G4ExtrudedSolid::DistanceToOut (const G4ThreeVector &p,
                                         const G4ThreeVector &v,
                                         const G4bool calcNorm,
                                               G4bool* validNorm,
                                               G4ThreeVector* n) const
{
  G4bool getnorm = calcNorm;
  if (getnorm) *validNorm = true;

  G4double z0 = fZSections[0].fZ;
  G4double z1 = fZSections[fNz-1].fZ;
  if ((p.z() <= z0 + kCarToleranceHalf) && v.z() < 0)
  {
    if (getnorm) n->set(0,0,-1);
    return 0;
  }
  if ((p.z() >= z1 - kCarToleranceHalf) && v.z() > 0)
  {
    if (getnorm) n->set(0,0,1);
    return 0;
  }

  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      // Intersection with Z planes
      //
      G4double dz = (z1 - z0)*0.5;
      G4double pz = p.z() - 0.5 * (z0 + z1);

      G4double vz = v.z();
      G4double tmax = (vz == 0) ? DBL_MAX : (std::copysign(dz,vz) - pz)/vz;
      G4int iside = (vz < 0) ? -4 : -2; // little trick: (-4+3)=-1, (-2+3)=+1

      // Intersection with lateral planes
      //
      std::size_t np = fPlanes.size();
      for (std::size_t i=0; i<np; ++i)
      {
        G4double cosa = fPlanes[i].a*v.x()+fPlanes[i].b*v.y();
        if (cosa > 0)
        {
          G4double dist = fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].d;
          if (dist >= -kCarToleranceHalf)
          {
            if (getnorm) n->set(fPlanes[i].a, fPlanes[i].b, fPlanes[i].c);
            return 0;
          }
          G4double tmp = -dist/cosa;
          if (tmax > tmp) { tmax = tmp; iside = (G4int)i; }
        }
      }

      // Set normal, if required, and return distance
      //
      if (getnorm)
      {
        if (iside < 0)
          { n->set(0, 0, iside + 3); } // (-4+3)=-1, (-2+3)=+1
        else
          { n->set(fPlanes[iside].a, fPlanes[iside].b, fPlanes[iside].c); }
      }
      return tmax;
    }
    case 2: // non-convex right prism
    {
    }
  }

  // Override the base class function to redefine validNorm
  // (the solid can be concave)

  G4double distOut =
    G4TessellatedSolid::DistanceToOut(p, v, calcNorm, validNorm, n);
  if (validNorm) { *validNorm = fIsConvex; }

  return distOut;
}

//_____________________________________________________________________________

G4double G4ExtrudedSolid::DistanceToOut(const G4ThreeVector& p) const
{
  switch (fSolidType)
  {
    case 1: // convex right prism
    {
      G4double dist = std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
      std::size_t np = fPlanes.size();
      for (std::size_t i=0; i<np; ++i)
      {
        G4double dd = fPlanes[i].a*p.x() + fPlanes[i].b*p.y() + fPlanes[i].d;
        if (dd > dist) dist = dd;
      }
      return (dist < 0) ? -dist : 0.;
    }
    case 2: // non-convex right prism
    {
      G4double distz = std::max(fZSections[0].fZ-p.z(),p.z()-fZSections[1].fZ);
      G4bool in = PointInPolygon(p);
      if (distz >= 0 || (!in)) return 0; // point is outside
      return std::min(-distz,std::sqrt(DistanceToPolygonSqr(p)));
    }
  }

  // General case: use tessellated solid
  return G4TessellatedSolid::DistanceToOut(p);
}

//_____________________________________________________________________________
// Get bounding box

void G4ExtrudedSolid::BoundingLimits(G4ThreeVector& pMin,
                                     G4ThreeVector& pMax) const
{
  G4double xmin0 = kInfinity, xmax0 = -kInfinity;
  G4double ymin0 = kInfinity, ymax0 = -kInfinity;

  for (G4int i=0; i<GetNofVertices(); ++i)
  {
    G4double x = fPolygon[i].x();
    if (x < xmin0) xmin0 = x;
    if (x > xmax0) xmax0 = x;
    G4double y = fPolygon[i].y();
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
    G4Exception("G4ExtrudedSolid::BoundingLimits()",
                "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//_____________________________________________________________________________
// Calculate extent under transform and specified limit

G4bool
G4ExtrudedSolid::CalculateExtent(const EAxis pAxis,
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
  G4TwoVectorList triangles;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  // triangulate the base polygon
  if (!G4GeomTools::TriangulatePolygon(fPolygon,triangles))
  {
    std::ostringstream message;
    message << "Triangulation of the base polygon has failed for solid: "
            << GetName() << " !"
            << "\nExtent has been calculated using boundary box";
    G4Exception("G4ExtrudedSolid::CalculateExtent()",
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
  G4int ntria = (G4int)triangles.size()/3;
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

//_____________________________________________________________________________

std::ostream& G4ExtrudedSolid::StreamInfo(std::ostream &os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid geometry type: " << fGeometryType  << G4endl;

  if ( fIsConvex) 
    { os << " Convex polygon; list of vertices:" << G4endl; }
  else  
    { os << " Concave polygon; list of vertices:" << G4endl; }
  
  for ( std::size_t i=0; i<fNv; ++i )
  {
    os << std::setw(5) << "#" << i 
       << "   vx = " << fPolygon[i].x()/mm << " mm" 
       << "   vy = " << fPolygon[i].y()/mm << " mm" << G4endl;
  }
  
  os << " Sections:" << G4endl;
  for ( std::size_t iz=0; iz<fNz; ++iz ) 
  {
    os << "   z = "   << fZSections[iz].fZ/mm          << " mm  "
       << "  x0= "    << fZSections[iz].fOffset.x()/mm << " mm  "
       << "  y0= "    << fZSections[iz].fOffset.y()/mm << " mm  " 
       << "  scale= " << fZSections[iz].fScale << G4endl;
  }     

/*
  // Triangles (for debugging)
  os << G4endl; 
  os << " Triangles:" << G4endl;
  os << " Triangle #   vertex1   vertex2   vertex3" << G4endl;

  G4int counter = 0;
  std::vector< std::vector<G4int> >::const_iterator it;
  for ( it = fTriangles.begin(); it != fTriangles.end(); it++ ) {
     std::vector<G4int> triangle = *it;
     os << std::setw(10) << counter++ 
        << std::setw(10) << triangle[0] << std::setw(10)  << triangle[1]
        << std::setw(10)  << triangle[2] 
        << G4endl;
  }          
*/
  os.precision(oldprc);

  return os;
}  

#endif
