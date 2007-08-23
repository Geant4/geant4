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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * and is subject to DEFCON 705 IPR conditions.                     *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.cc,v 1.10 2007-08-23 14:49:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4TessellatedSolid.cc
//
// Date:                15/06/2005
// Author:              P R Truscott
// Organisation:        QinetiQ Ltd, UK
// Customer:            UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:            C/MAT/N03517
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 22 November 2005, F Lei
//  - Changed ::DescribeYourselfTo(), line 464
//  - added GetPolyHedron()
// 
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "G4TessellatedSolid.hh"
#include "G4PolyhedronArbitrary.hh"
#include "globals.hh"
#include "Randomize.hh"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no facets.
//
G4TessellatedSolid::G4TessellatedSolid ()
  : G4VSolid("dummy"), fpPolyhedron(0), cubicVolume(0.), surfaceArea(0.)
{
  dirTolerance = 1.0E-14;
  
  geometryType = "G4TessellatedSolid";
  facets.clear();
  solidClosed  = false;
  
  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;
}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor. Simple define name and geometry type - no facets
// to detine.
//
G4TessellatedSolid::G4TessellatedSolid (const G4String &name)
  : G4VSolid(name), fpPolyhedron(0), cubicVolume(0.), surfaceArea(0.)
{
  dirTolerance = 1.0E-14;
  
  geometryType = "G4TessellatedSolid";
  facets.clear();
  solidClosed  = false;
  
  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;
}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4TessellatedSolid::G4TessellatedSolid( __void__& a )
  : G4VSolid(a), fpPolyhedron(0), facets(0),
    geometryType("G4TessellatedSolid"), cubicVolume(0.), surfaceArea(0.),
    vertexList(), xMinExtent(0.), xMaxExtent(0.),
    yMinExtent(0.), yMaxExtent(0.), zMinExtent(0.), zMaxExtent(0.),
    solidClosed(false)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor.
//
G4TessellatedSolid::~G4TessellatedSolid ()
{
  DeleteObjects ();
}

///////////////////////////////////////////////////////////////////////////////
//
// Define copy constructor.
//
G4TessellatedSolid::G4TessellatedSolid (const G4TessellatedSolid &s)
  : G4VSolid(s)
{
  if (&s == this) { return; }

  CopyObjects (s);
}

///////////////////////////////////////////////////////////////////////////////
//
// Define assignment operator.
//
const G4TessellatedSolid &
G4TessellatedSolid::operator= (const G4TessellatedSolid &s)
{
  if (&s == this) { return *this; }
  
  DeleteObjects ();
  CopyObjects (s);
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DeleteObjects ()
{
  for (std::vector<G4VFacet *>::iterator f=facets.begin(); f!=facets.end(); f++)
  {
    delete *f;
  }
  facets.clear();
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::CopyObjects (const G4TessellatedSolid &s)
{
  size_t n = s.GetNumberOfFacets();
  for (size_t i=0; i<n; n++)
  {
    G4VFacet *facetClone = (s.GetFacet(i))->GetClone();
    AddFacet(facetClone);
  }
  solidClosed = s.GetSolidClosed();

//  cubicVolume = s.GetCubicVolume();  
}

///////////////////////////////////////////////////////////////////////////////
//
// Add a facet to the facet list.  Note that you can add, but you cannot
// delete.
//
G4bool G4TessellatedSolid::AddFacet (G4VFacet *aFacet)
{
  // Add the facet to the vector.

  if (solidClosed)
  {
    G4Exception("G4TessellatedSolid::AddFacet()", "InvalidSetup",
                JustWarning, "Attempt to add facets when solid is closed.");
    return false;
  }
  else if (aFacet->IsDefined())
  {
    if (facets.size() == 0)
    {
      facets.push_back(aFacet);
    }
    else
    {
      G4bool found = false;
      FacetI it    = facets.begin();
      do
      {
        found = (**it == *aFacet);
      } while (!found && ++it!=facets.end());
    
      if (found)
      {
        delete *it;
        facets.erase(it);
      }
      else
      {
        facets.push_back(aFacet);
      }
    }
    
    return true;
  }
  else
  {
    G4Exception("G4TessellatedSolid::AddFacet()", "InvalidSetup",
                JustWarning, "Attempt to add facet not properly defined.");
    G4cerr << "Facet attributes:" << G4endl;
    aFacet->StreamInfo(G4cerr);
    G4cerr << G4endl;
    
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::SetSolidClosed (const G4bool t)
{
  if (t)
  {
    vertexList.clear();
    for (FacetCI it=facets.begin(); it!=facets.end(); it++)
    {
      size_t m = vertexList.size();
      G4ThreeVector p(0.0,0.0,0.0);
      for (size_t i=0; i<(*it)->GetNumberOfVertices(); i++)
      {
        p            = (*it)->GetVertex(i);
        G4bool found = false;
        size_t j     = 0;
        while (j < m && !found)
        {
          G4ThreeVector q = vertexList[j];
          found = (q-p).mag() < 0.5*kCarTolerance;
          if (!found) j++;
        }
    
        if (!found)
        {
          vertexList.push_back(p);
          (*it)->SetVertexIndex(i,vertexList.size()-1);
        }
        else
        {
          (*it)->SetVertexIndex(i,j);
        }
      }
    }
    //
    // Now update the maximum x, y and z limits of the volume.
    //
    for (size_t i=0; i<vertexList.size(); i++)
    {
      G4ThreeVector p = vertexList[i];
      G4double x      = p.x();
      G4double y      = p.y();
      G4double z      = p.z();
    
      if (i > 0)
      {
        if (x > xMaxExtent) xMaxExtent = x;
        if (x < xMinExtent) xMinExtent = x;
        if (y > yMaxExtent) yMaxExtent = y;
        if (y < yMinExtent) yMinExtent = y;
        if (z > zMaxExtent) zMaxExtent = z;
        if (z < zMinExtent) zMinExtent = z;
      }
      else
      {
        xMaxExtent = x;
        xMinExtent = x;
        yMaxExtent = y;
        yMinExtent = y;
        zMaxExtent = z;
        zMinExtent = z;    
      }
    }
//
//
// Compute extremeFacets, i.e. find those facets that have surface
// planes that bound the volume.
// Note that this is going to reject concaved surfaces as being extreme.  Also
// note that if the vertex is on the facet, displacement is zero, so IsInside
// returns true.  So will this work??  Need non-equality
// "G4bool inside = displacement < 0.0;"
// or
// "G4bool inside = displacement <= -0.5*kCarTolerance" 
// (Notes from PT 13/08/2007).
//
    for (FacetCI it=facets.begin(); it!=facets.end(); it++)
    {
      G4bool isExtreme = true;
      for (size_t i=0; i<vertexList.size(); i++)
      {
        if (!(*it)->IsInside(vertexList[i]))
        {
          isExtreme = false;
          break;
        }
      }
      if (isExtreme)
        extremeFacets.insert(*it);
    }
    solidClosed = true;
  }
  else
  {
    solidClosed = false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// GetSolidClosed
//
// Used to determine whether the solid is closed to adding further facets.
//
G4bool G4TessellatedSolid::GetSolidClosed () const
  {return solidClosed;}

///////////////////////////////////////////////////////////////////////////////
//
// operator+=
//
// This operator allows the user to add two tessellated solids together, so
// that the solid on the left then includes all of the facets in the solid
// on the right.  Note that copies of the facets are generated, rather than
// using the original facet set of the solid on the right.
//
const G4TessellatedSolid &G4TessellatedSolid::operator+=
  (const G4TessellatedSolid &right)
{
  for (size_t i=0; i<right.GetNumberOfFacets(); i++)
    AddFacet(right.GetFacet(i)->GetClone());
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetFacet
//
// Access pointer to facet in solid, indexed by integer i.
//
G4VFacet *G4TessellatedSolid::GetFacet (size_t i) const
{
  return facets[i];
}

///////////////////////////////////////////////////////////////////////////////
//
// GetNumberOfFacets
//
size_t G4TessellatedSolid::GetNumberOfFacets () const
{
  return facets.size();
}

///////////////////////////////////////////////////////////////////////////////
//
// EInside G4TessellatedSolid::Inside (const G4ThreeVector &p) const
//
// This method must return:
//    * kOutside if the point at offset p is outside the shape
//      boundaries plus kCarTolerance/2,
//    * kSurface if the point is <= kCarTolerance/2 from a surface, or
//    * kInside otherwise.
//
EInside G4TessellatedSolid::Inside (const G4ThreeVector &p) const
{
  if ( p.x() < xMinExtent - kCarTolerance ||
       p.x() > xMaxExtent + kCarTolerance ||
       p.y() < yMinExtent - kCarTolerance ||
       p.y() > yMaxExtent + kCarTolerance ||
       p.z() < zMinExtent - kCarTolerance ||
       p.z() > zMaxExtent + kCarTolerance )
  {
    return kOutside;
  }  

  G4double minDist = kInfinity;
//
//
// Check if we are close to a surface
//
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    G4double dist = (*f)->Distance(p,minDist);
    if (dist < minDist) minDist = dist;
    if (dist <= 0.5*kCarTolerance)
    {
      return kSurface;
    }
  }
//
//
// The following is something of an adaptation of the method implemented by
// Rickard Holmberg augmented with information from Schneider & Eberly,
// "Geometric Tools for Computer Graphics," pp700-701, 2003.  In essence, we're
// trying to determine whether we're inside the volume by projecting a few rays
// and determining if the first surface crossed is has a normal vector between
// 0 to pi/2 (out-going) or pi/2 to pi (in-going).  We should also avoid rays
// which are nearly within the plane of the tessellated surface, and therefore
// produce rays randomly.  For the moment, this is a bit over-engineered
// (belt-braces-and-ducttape).
//
  G4int nTry                = 7;
  G4double distO            = 0.0;
  G4double distI            = 0.0;
  G4double distFromSurfaceO = 0.0;
  G4double distFromSurfaceI = 0.0;
  G4ThreeVector normalO(0.0,0.0,0.0);
  G4ThreeVector normalI(0.0,0.0,0.0);
  G4bool crossingO          = false;
  G4bool crossingI          = false;
  EInside location          = kOutside;
  EInside locationprime     = kOutside;

  for (G4int i=0; i<nTry; i++)
  {
    G4double distOut = kInfinity;
    G4double distIn  = kInfinity;
    G4bool nearParallel = false;
    do
    {
      distOut          = kInfinity;
      distIn           = kInfinity;
      G4ThreeVector v  = G4ThreeVector(G4UniformRand()-0.5,
        G4UniformRand()-0.5, G4UniformRand()-0.5).unit();
      FacetCI f = facets.begin();
      do
      {
        crossingO =  ((*f)->Intersect(p,v,true,distO,distFromSurfaceO,normalO));
        crossingI =  ((*f)->Intersect(p,v,false,distI,distFromSurfaceI,normalI));
        if (crossingO || crossingI)
        {
          nearParallel = crossingO && std::abs(normalO.dot(v))<dirTolerance ||
                         crossingI && std::abs(normalI.dot(v))<dirTolerance;
          if (!nearParallel)
          {
            if (crossingO && distO > 0.0 && distO < distOut) distOut = distO;
            if (crossingI && distI > 0.0 && distI < distIn)  distIn  = distI;
          }
        }
      } while (!nearParallel && ++f!=facets.end());
    } while (nearParallel);
    if (distIn == kInfinity && distOut == kInfinity)
      locationprime = kOutside;
    else if (distIn <= distOut - kCarTolerance*0.5)
      locationprime = kOutside;
    else if (distOut <= distIn - kCarTolerance*0.5)
      locationprime = kInside;

    if (i == 0) location = locationprime;
    else if (locationprime != location)
    {
      G4Exception("G4TessellatedSolid::Inside()()",
                  "UnknownInsideOutside", FatalException,
                  "Cannot determine whether point is inside or outside volume !" );
    }
  }

  return location;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4ThreeVector G4TessellatedSolid::SurfaceNormal (const G4ThreeVector &p) const
//
// Return the outwards pointing unit normal of the shape for the
// surface closest to the point at offset p.

G4ThreeVector G4TessellatedSolid::SurfaceNormal (const G4ThreeVector &p) const
{
  FacetCI minFacet;
  G4double minDist   = kInfinity;
  G4double dist      = 0.0;
  G4ThreeVector normal;
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist);
    if (dist < minDist)
    {
      minDist  = dist;
      minFacet = f;
    }
  }
  
  if (minDist != kInfinity)
  {
     normal = (*minFacet)->GetSurfaceNormal();
  }
  else
  {
#ifdef G4VERBOSE
    G4cout << "WARNING - G4TessellatedSolid::SurfaceNormal(p)" << G4endl
           << "          No facets found for point: " << p << " !" << G4endl
           << "          Returning approximated value for normal." << G4endl;
    G4Exception("G4TessellatedSolid::SurfaceNormal(p)", "Notification",
                JustWarning, "Point p is not on surface !?" );
#endif
    normal = (p.z()>0 ? G4ThreeVector(0,0,1) : G4ThreeVector(0,0,-1));
  }

  return normal;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v)
//
// Return the distance along the normalised vector v to the shape,
// from the point at offset p. If there is no intersection, return
// kInfinity. The first intersection resulting from ‘leaving’ a
// surface/volume is discarded. Hence, this is tolerant of points on
// surface of shape.

G4double G4TessellatedSolid::DistanceToIn (const G4ThreeVector &p,
  const G4ThreeVector &v) const
{
  G4double minDist         = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal(0.0,0.0,0.0);
  
#if G4SPECSDEBUG
  if ( Inside(p) == kInside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     //     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout << "DistanceToOut(p) == " << DistanceToOut(p) << G4endl;
     G4Exception("G4TriangularFacet::DistanceToIn(p,v)", "Notification", JustWarning, 
                 "Point p is already inside!?" );
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    if ((*f)->Intersect(p,v,false,dist,distFromSurface,normal))
    {
      if (distFromSurface > 0.5*kCarTolerance && dist >= 0.0 &&
        dist < minDist) minDist  = dist;
    }
  }

  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToIn(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an outside point p. The
// distance can be an underestimate.

G4double G4TessellatedSolid::DistanceToIn (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double dist    = 0.0;
  
#if G4SPECSDEBUG
  if ( Inside(p) == kInside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     //     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout << "DistanceToOut(p) == " << DistanceToOut(p) << G4endl;
     G4Exception("G4TriangularFacet::DistanceToIn(p)", "Notification", JustWarning, 
                 "Point p is already inside!?" );
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist,false);
    if (dist < minDist) minDist  = dist;
  }
  
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
//                        const G4bool calcNorm=false,
//                        G4bool *validNorm=0, G4ThreeVector *n=0);
//
// Return distance along the normalised vector v to the shape, from a
// point at an offset p inside or on the surface of the
// shape. Intersections with surfaces, when the point is not greater
// than kCarTolerance/2 from a surface, must be ignored.
//     If calcNorm is true, then it must also set validNorm to either
//     * true, if the solid lies entirely behind or on the exiting
//        surface. Then it must set n to the outwards normal vector
//        (the Magnitude of the vector is not defined).
//     * false, if the solid does not lie entirely behind or on the
//       exiting surface.
// If calcNorm is false, then validNorm and n are unused.

G4double G4TessellatedSolid::DistanceToOut (const G4ThreeVector &p,
                    const G4ThreeVector &v, const G4bool calcNorm,
                          G4bool *validNorm, G4ThreeVector *n) const
{
  G4double minDist         = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal(0.0,0.0,0.0);
  G4ThreeVector minNormal(0.0,0.0,0.0);
  
#if G4SPECSDEBUG
  if ( Inside(p) == kOutside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     //     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout << "DistanceToIn(p) == " << DistanceToIn(p) << G4endl;
     G4Exception("G4TriangularFacet::DistanceToOut(p)", "Notification", JustWarning, 
                 "Point p is already outside !?" );
  }
#endif

  G4bool isExtreme = false;
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    if ((*f)->Intersect(p,v,true,dist,distFromSurface,normal))
     {
      if (distFromSurface > 0.0 && distFromSurface <= 0.5*kCarTolerance &&
          (*f)->Distance(p,kCarTolerance) <= 0.5*kCarTolerance)
      {
        // We are on a surface. Return zero.
        *validNorm = extremeFacets.count(*f);
        *n         = SurfaceNormal(p);
        return 0.0;
      }
      if (dist >= 0.0 && dist < minDist)
      {
        minDist   = dist;
        minNormal = normal;
        isExtreme = extremeFacets.count(*f);
      }
    }
  }
  
  if (minDist < kInfinity)
  {
    if (calcNorm)
    {
      *validNorm = isExtreme;
      *n         = minNormal;
    }
    return minDist;
  }
  else
  {
    // No intersection found
    if (calcNorm)
    {
      *validNorm = false;
      *n         = SurfaceNormal(p);
    }
    return 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// G4double DistanceToOut(const G4ThreeVector& p)
//
// Calculate distance to nearest surface of shape from an inside
// point. The distance can be an underestimate.

G4double G4TessellatedSolid::DistanceToOut (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double dist    = 0.0;
  
#if G4SPECSDEBUG
  if ( Inside(p) == kOutside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     //     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout << "DistanceToIn(p) == " << DistanceToIn(p) << G4endl;
     G4Exception("G4TriangularFacet::DistanceToOut(p)", "Notification", JustWarning, 
                 "Point p is already outside !?" );
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist,true);
    if (dist < minDist) minDist  = dist;
  }
  
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
// G4GeometryType GetEntityType() const;
//
// Provide identification of the class of an object (required for
// persistency and STEP interface).
//
G4GeometryType G4TessellatedSolid::GetEntityType () const
{
  return geometryType;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

///////////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//                                                                                
//void G4TessellatedSolid::ComputeDimensions (G4VPVParameterisation* p,
//  const G4int n, const G4VPhysicalVolume* pRep) const
//{
//  G4VSolid *ptr = 0;
//  ptr           = *this;
//  p->ComputeDimensions(ptr,n,pRep);
//}

///////////////////////////////////////////////////////////////////////////////
//
std::ostream &G4TessellatedSolid::StreamInfo(std::ostream &os) const
{
  os << G4endl;
  os << "Geometry Type    = " << geometryType  << G4endl;
  os << "Number of facets = " << facets.size() << G4endl;
  
  for (FacetCI f = facets.begin(); f != facets.end(); f++)
  {
    os << "FACET #          = " << f-facets.begin()+1 << G4endl;
    (*f)->StreamInfo(os);
  }
  os <<G4endl;
  
  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
G4Polyhedron *G4TessellatedSolid::CreatePolyhedron () const
{
  size_t nVertices = vertexList.size();
  size_t nFacets   = facets.size();
  G4PolyhedronArbitrary *polyhedron =
    new G4PolyhedronArbitrary (nVertices, nFacets);
  for (G4ThreeVectorList::const_iterator v = vertexList.begin();
        v!=vertexList.end(); v++) polyhedron->AddVertex(*v);
    
  for (FacetCI f=facets.begin(); f != facets.end(); f++)
  {
    size_t v[4];
    for (size_t j=0; j<4; j++)
    {
      size_t i = (*f)->GetVertexIndex(j);
      if (i == 999999999) v[j] = 0;
      else                v[j] = i+1;
    }
    if ((*f)->GetEntityType() == "G4RectangularFacet")
    {
      size_t i = v[3];
      v[3]     = v[2];
      v[2]     = i;
    }
    polyhedron->AddFacet(v[0],v[1],v[2],v[3]);
  }
  
  return (G4Polyhedron*) polyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
G4NURBS *G4TessellatedSolid::CreateNURBS () const
{
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4TessellatedSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4TessellatedSolid::CalculateExtent(const EAxis pAxis,
         const G4VoxelLimits& pVoxelLimit, const G4AffineTransform& pTransform,
               G4double& pMin, G4double& pMax) const
{
  if (!pTransform.IsRotated())
  {
    G4double xoffset,xMin,xMax;
    G4double yoffset,yMin,yMax;
    G4double zoffset,zMin,zMax;
                                                                                
    xoffset = pTransform.NetTranslation().x();
    xMin    = xoffset + xMinExtent;
    xMax    = xoffset + xMaxExtent;
                                                                                
    if (pVoxelLimit.IsXLimited())
    {
      if ( xMin > pVoxelLimit.GetMaxXExtent()+kCarTolerance ||
           xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance    ) return false ;
      else
      {
        if (xMin < pVoxelLimit.GetMinXExtent())
        {
          xMin = pVoxelLimit.GetMinXExtent() ;
        }
        if (xMax > pVoxelLimit.GetMaxXExtent())
        {
          xMax = pVoxelLimit.GetMaxXExtent() ;
        }
      }
    }
 
    yoffset = pTransform.NetTranslation().y();
    yMin    = yoffset + yMinExtent;
    yMax    = yoffset + yMaxExtent;
                                                                                
    if (pVoxelLimit.IsYLimited())
    {
      if ( yMin > pVoxelLimit.GetMaxYExtent()+kCarTolerance ||
           yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance    ) return false ;
      else
      {
        if (yMin < pVoxelLimit.GetMinYExtent())
        {
          yMin = pVoxelLimit.GetMinYExtent() ;
        }
        if (yMax > pVoxelLimit.GetMaxYExtent())
        {
          yMax = pVoxelLimit.GetMaxYExtent() ;
        }
      }
    }

    zoffset = pTransform.NetTranslation().z();
    zMin    = zoffset + zMinExtent;
    zMax    = zoffset + zMaxExtent;
                                                                                
    if (pVoxelLimit.IsZLimited())
    {
      if ( zMin > pVoxelLimit.GetMaxZExtent()+kCarTolerance ||
           zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance    ) return false ;
      else
      {
        if (zMin < pVoxelLimit.GetMinZExtent())
        {
          zMin = pVoxelLimit.GetMinZExtent() ;
        }
        if (zMax > pVoxelLimit.GetMaxZExtent())
        {
          zMax = pVoxelLimit.GetMaxZExtent() ;
        }
      }
    }

    switch (pAxis)
    {
      case kXAxis:
        pMin = xMin ;
        pMax = xMax ;
        break ;
      case kYAxis:
        pMin=yMin;
        pMax=yMax;
        break;
      case kZAxis:
        pMin=zMin;
        pMax=zMax;
        break;
      default:
        break;
    }
    pMin -= kCarTolerance ;
    pMax += kCarTolerance ;
                                                                                
    return true;
  }
  else
  {
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinXExtent () const
  {return xMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxXExtent () const
  {return xMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinYExtent () const
  {return yMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxYExtent () const
  {return yMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMinZExtent () const
  {return zMinExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetMaxZExtent () const
  {return zMaxExtent;}

///////////////////////////////////////////////////////////////////////////////
//
G4VisExtent G4TessellatedSolid::GetExtent () const
{
  return G4VisExtent (xMinExtent, xMaxExtent, yMinExtent, yMaxExtent,
    zMinExtent, zMaxExtent);
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetCubicVolume ()
{
  if(cubicVolume != 0.) {;}
  else   { cubicVolume = G4VSolid::GetCubicVolume(); }
  return cubicVolume;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::GetSurfaceArea ()
{
  if(surfaceArea != 0.) { return surfaceArea; }

  for (FacetI f=facets.begin(); f!=facets.end(); f++)
  {
    surfaceArea += (*f)->GetArea();
  }
  return surfaceArea;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::GetPointOnSurface() const
{
  // Select randomly a facet and return a random point on it
 
  G4int i = CLHEP::RandFlat::shootInt(facets.size());
  return facets[i]->GetPointOnFace();
}
