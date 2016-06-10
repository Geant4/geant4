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
//
// $Id: G4ExtrudedSolid.cc 95301 2016-02-04 11:25:33Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ExtrudedSolid.cc
//
// Author: Ivana Hrivnacova, IPN Orsay
// --------------------------------------------------------------------

#include "G4ExtrudedSolid.hh"

#if !defined(G4GEOM_USE_UEXTRUDEDSOLID)

#include <set>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VFacet.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( const G4String& pName,
                                        std::vector<G4TwoVector> polygon,
                                        std::vector<ZSection> zsections)
  : G4TessellatedSolid(pName),
    fNv(polygon.size()),
    fNz(zsections.size()),
    fPolygon(),
    fZSections(),
    fTriangles(),
    fIsConvex(false),
    fGeometryType("G4ExtrudedSolid")
    
{
  // General constructor 

  // First check input parameters

  if ( fNv < 3 )
  {
    std::ostringstream message;
    message << "Number of polygon vertices < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  if ( fNz < 2 )
  {
    std::ostringstream message;
    message << "Number of z-sides < 2 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  for ( G4int i=0; i<fNz-1; ++i ) 
  {
    if ( zsections[i].fZ > zsections[i+1].fZ ) 
    {
      std::ostringstream message;
      message << "Z-sections have to be ordered by z value (z0 < z1 < z2...) - "
              << pName;
      G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                  FatalErrorInArgument, message);
    }
    if ( std::fabs( zsections[i+1].fZ - zsections[i].fZ ) < kCarTolerance * 0.5 ) 
    {
      std::ostringstream message;
      message << "Z-sections with the same z position are not supported - "
              << pName;
      G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0001",
                  FatalException, message);
    }
  }  
  
  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)
  //
  G4double area = 0.;
  for ( G4int i=0; i<fNv; ++i ) {
    G4int j = i+1;
    if ( j == fNv ) j = 0;
    area += 0.5 * ( polygon[i].x()*polygon[j].y() - polygon[j].x()*polygon[i].y());
  }
 
  // Copy polygon
  //
  if  ( area < 0. ) {   
    // Polygon vertices are defined clockwise, we just copy the polygon       
    for ( G4int i=0; i<fNv; ++i ) { fPolygon.push_back(polygon[i]); }
  }
  else {
    // Polygon vertices are defined anti-clockwise, we revert them
    //G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
    //            JustWarning, 
    //            "Polygon vertices defined anti-clockwise, reverting polygon");      
    for ( G4int i=0; i<fNv; ++i ) { fPolygon.push_back(polygon[fNv-i-1]); }
  }
    
  
  // Copy z-sections
  //
  for ( G4int i=0; i<fNz; ++i ) { fZSections.push_back(zsections[i]); }
    

  G4bool result = MakeFacets();
  if (!result)
  {   
    std::ostringstream message;
    message << "Making facets failed - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0003",
                FatalException, message);
  }
  fIsConvex = IsConvex();

  
  ComputeProjectionParameters();
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( const G4String& pName,
                                        std::vector<G4TwoVector> polygon,       
                                        G4double dz,
                                        G4TwoVector off1, G4double scale1,
                                        G4TwoVector off2, G4double scale2 )
  : G4TessellatedSolid(pName),
    fNv(polygon.size()),
    fNz(2),
    fPolygon(),
    fZSections(),
    fTriangles(),
    fIsConvex(false),
    fGeometryType("G4ExtrudedSolid")
    
{
  // Special constructor for solid with 2 z-sections

  // First check input parameters
  //
  if ( fNv < 3 )
  {
    std::ostringstream message;
    message << "Number of polygon vertices < 3 - " << pName;
    G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
     
  // Check if polygon vertices are defined clockwise
  // (the area is positive if polygon vertices are defined anti-clockwise)
  
  G4double area = 0.;
  for ( G4int i=0; i<fNv; ++i )
  {
    G4int j = i+1;
    if ( j == fNv )  { j = 0; }
    area += 0.5 * ( polygon[i].x()*polygon[j].y()
                  - polygon[j].x()*polygon[i].y());
  }
 
  // Copy polygon
  //
  if  ( area < 0. )
  {   
    // Polygon vertices are defined clockwise, we just copy the polygon       
    for ( G4int i=0; i<fNv; ++i ) { fPolygon.push_back(polygon[i]); }
  }
  else
  {
    // Polygon vertices are defined anti-clockwise, we revert them
    //G4Exception("G4ExtrudedSolid::G4ExtrudedSolid()", "GeomSolids1001",
    //            JustWarning, 
    //            "Polygon vertices defined anti-clockwise, reverting polygon");      
    for ( G4int i=0; i<fNv; ++i ) { fPolygon.push_back(polygon[fNv-i-1]); }
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
  fIsConvex = IsConvex();

  ComputeProjectionParameters();
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid( __void__& a )
  : G4TessellatedSolid(a), fNv(0), fNz(0), fPolygon(), fZSections(),
    fTriangles(), fIsConvex(false), fGeometryType("G4ExtrudedSolid")
{
  // Fake default constructor - sets only member data and allocates memory
  //                            for usage restricted to object persistency.
}

//_____________________________________________________________________________

G4ExtrudedSolid::G4ExtrudedSolid(const G4ExtrudedSolid& rhs)
  : G4TessellatedSolid(rhs), fNv(rhs.fNv), fNz(rhs.fNz),
    fPolygon(rhs.fPolygon), fZSections(rhs.fZSections),
    fTriangles(rhs.fTriangles), fIsConvex(rhs.fIsConvex),
    fGeometryType(rhs.fGeometryType), fKScales(rhs.fKScales),
    fScale0s(rhs.fScale0s), fKOffsets(rhs.fKOffsets), fOffset0s(rhs.fOffset0s)
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
   fGeometryType = rhs.fGeometryType; fKScales = rhs.fKScales;
   fScale0s = rhs.fScale0s; fKOffsets = rhs.fKOffsets;
   fOffset0s = rhs.fOffset0s;

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

  for ( G4int iz=0; iz<fNz-1; ++iz) 
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
  G4int iz = 0;
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

G4bool G4ExtrudedSolid::IsSameLine(G4TwoVector p,  
                                   G4TwoVector l1, G4TwoVector l2) const
{
  // Return true if p is on the line through l1, l2 

  if ( l1.x() == l2.x() )
  {
    return std::fabs(p.x() - l1.x()) < kCarTolerance * 0.5; 
  }
   G4double  slope= ((l2.y() - l1.y())/(l2.x() - l1.x())); 
   G4double predy= l1.y() +  slope *(p.x() - l1.x()); 
   G4double dy= p.y() - predy; 

   // Calculate perpendicular distance
   //
   // G4double perpD= std::fabs(dy) / std::sqrt( 1 + slope * slope ); 
   // G4bool   simpleComp= (perpD<0.5*kCarTolerance); 

   // Check perpendicular distance vs tolerance 'directly'
   //
   const G4double tol= 0.5 * kCarTolerance ; 
   G4bool    squareComp=  (dy*dy < (1+slope*slope) * tol * tol); 

   // return  simpleComp; 
   return squareComp;
}                    

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsSameLineSegment(G4TwoVector p,  
                                   G4TwoVector l1, G4TwoVector l2) const
{
  // Return true if p is on the line through l1, l2 and lies between
  // l1 and l2 

  if ( p.x() < std::min(l1.x(), l2.x()) - kCarTolerance * 0.5 || 
       p.x() > std::max(l1.x(), l2.x()) + kCarTolerance * 0.5 ||
       p.y() < std::min(l1.y(), l2.y()) - kCarTolerance * 0.5 || 
       p.y() > std::max(l1.y(), l2.y()) + kCarTolerance * 0.5 )
  {
    return false;
  }

  return IsSameLine(p, l1, l2);
}

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsSameSide(G4TwoVector p1, G4TwoVector p2, 
                                   G4TwoVector l1, G4TwoVector l2) const
{
  // Return true if p1 and p2 are on the same side of the line through l1, l2 

  return   ( (p1.x() - l1.x()) * (l2.y() - l1.y())
         - (l2.x() - l1.x()) * (p1.y() - l1.y()) )
         * ( (p2.x() - l1.x()) * (l2.y() - l1.y())
         - (l2.x() - l1.x()) * (p2.y() - l1.y()) ) > 0;
}       

//_____________________________________________________________________________

G4bool G4ExtrudedSolid::IsPointInside(G4TwoVector a, G4TwoVector b,
                                      G4TwoVector c, G4TwoVector p) const
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
G4ExtrudedSolid::GetAngle(G4TwoVector po, G4TwoVector pa, G4TwoVector pb) const
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
    // vertices ardered clock wise has to be reordered

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
  vertices.push_back(GetVertex(fNz-1, ind1));
  vertices.push_back(GetVertex(fNz-1, ind2));
  vertices.push_back(GetVertex(fNz-1, ind3));
  
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

  // Fill one more vector
  //
  std::vector< Vertex > verticesToBeDone;
  for ( G4int i=0; i<fNv; ++i )
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

    G4int counter = 0;
    while ( angle >= pi )    // Loop checking, 13.08.2015, G.Cosmo
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
      
      counter++;
      
      if ( counter > fNv) {
        G4Exception("G4ExtrudedSolid::AddGeneralPolygonFacets",
                    "GeomSolids0003", FatalException,
                    "Triangularisation has failed.");
        break;
      }  
    }

    G4bool good = true;
    std::vector< Vertex >::iterator it;
    for ( it=verticesToBeDone.begin(); it != verticesToBeDone.end(); ++it )
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

    good = AddFacet( new G4TriangularFacet( GetVertex(fNz-1, 2), GetVertex(fNz-1, 1),
                                            GetVertex(fNz-1, 0), ABSOLUTE) );
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

    good = AddFacet( new G4QuadrangularFacet( GetVertex(fNz-1, 3), GetVertex(fNz-1, 2), 
                                              GetVertex(fNz-1, 1), GetVertex(fNz-1, 0),
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
  for ( G4int iz = 0; iz < fNz-1; ++iz ) 
  {
    for ( G4int i = 0; i < fNv; ++i )
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

G4bool G4ExtrudedSolid::IsConvex() const
{
  // Get polygon convexity (polygon is convex if all vertex angles are < pi )

  for ( G4int i=0; i< fNv; ++i )
  {
    G4int j = ( i + 1 ) % fNv;
    G4int k = ( i + 2 ) % fNv;
    G4TwoVector v1 = fPolygon[i]-fPolygon[j];
    G4TwoVector v2 = fPolygon[k]-fPolygon[j];
    G4double dphi = v2.phi() - v1.phi();
    if ( dphi < 0. )  { dphi += 2.*pi; }
    
    if ( dphi >= pi ) { return false; }
  }
  
  return true;
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

EInside G4ExtrudedSolid::Inside (const G4ThreeVector &p) const
{
  // Override the base class function  as it fails in case of concave polygon.
  // Project the point in the original polygon scale and check if it is inside
  // for each triangle.

  // Check first if outside extent
  //
  if ( p.x() < GetMinXExtent() - kCarTolerance * 0.5 ||
       p.x() > GetMaxXExtent() + kCarTolerance * 0.5 ||
       p.y() < GetMinYExtent() - kCarTolerance * 0.5 ||
       p.y() > GetMaxYExtent() + kCarTolerance * 0.5 ||
       p.z() < GetMinZExtent() - kCarTolerance * 0.5 ||
       p.z() > GetMaxZExtent() + kCarTolerance * 0.5 )
  {
    // G4cout << "G4ExtrudedSolid::Outside extent: " << p << G4endl;
    return kOutside;
  }  

  // Project point p(z) to the polygon scale p0
  //
  G4TwoVector pscaled = ProjectPoint(p);
  
  // Check if on surface of polygon
  //
  for ( G4int i=0; i<fNv; ++i )
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
  std::vector< std::vector<G4int> >::const_iterator it = fTriangles.begin();
  G4bool inside = false;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    if ( IsPointInside(fPolygon[(*it)[0]], fPolygon[(*it)[1]],
                       fPolygon[(*it)[2]], pscaled) )  { inside = true; }
    ++it;
  } while ( (inside == false) && (it != fTriangles.end()) );
  
  if ( inside )
  {
    // Check if on surface of z sides
    //
    if ( std::fabs( p.z() - fZSections[0].fZ ) < kCarTolerance * 0.5 ||
         std::fabs( p.z() - fZSections[fNz-1].fZ ) < kCarTolerance * 0.5 )
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

G4double G4ExtrudedSolid::DistanceToOut (const G4ThreeVector &p,
                                         const G4ThreeVector &v,
                                         const G4bool calcNorm,
                                               G4bool *validNorm,
                                               G4ThreeVector *n) const
{
  // Override the base class function to redefine validNorm
  // (the solid can be concave) 

  G4double distOut =
    G4TessellatedSolid::DistanceToOut(p, v, calcNorm, validNorm, n);
  if (validNorm) { *validNorm = fIsConvex; }

  return distOut;
}


//_____________________________________________________________________________

G4double G4ExtrudedSolid::DistanceToOut (const G4ThreeVector &p) const
{
  // Override the overloaded base class function

  return G4TessellatedSolid::DistanceToOut(p);
}

//_____________________________________________________________________________

std::ostream& G4ExtrudedSolid::StreamInfo(std::ostream &os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid geometry type: " << fGeometryType  << G4endl;

  if ( fIsConvex) 
    { os << " Convex polygon; list of vertices:" << G4endl; }
  else  
    { os << " Concave polygon; list of vertices:" << G4endl; }
  
  for ( G4int i=0; i<fNv; ++i )
  {
    os << std::setw(5) << "#" << i 
       << "   vx = " << fPolygon[i].x()/mm << " mm" 
       << "   vy = " << fPolygon[i].y()/mm << " mm" << G4endl;
  }
  
  os << " Sections:" << G4endl;
  for ( G4int iz=0; iz<fNz; ++iz ) 
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
        << std::setw(10) << triangle[0] << std::setw(10)  << triangle[1]  << std::setw(10)  << triangle[2] 
        << G4endl;
  }          
*/
  os.precision(oldprc);

  return os;
}  

#endif
