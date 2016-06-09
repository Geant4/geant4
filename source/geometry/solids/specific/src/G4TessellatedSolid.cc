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
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.cc,v 1.27 2010-11-02 11:29:07 gcosmo Exp $
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
// 22 August 2011,    I Hrivnacova, Orsay, fix in DistanceToOut(p) and
//                    DistanceToIn(p) to exactly compute distance from facet
//                    avoiding use of 'outgoing' flag shortcut variant.
//
// 04 August 2011,    T Nikitina, CERN, added SetReferences() to
//                    CreatePolyhedron() for Visualization of Boolean Operations  
//
// 12 April 2010,     P R Truscott, QinetiQ, bug fixes to treat optical
//                    photon transport, in particular internal reflection
//                    at surface.
//
// 14 November 2007,  P R Truscott, QinetiQ & Stan Seibert, U Texas
//                    Bug fixes to CalculateExtent
//
// 17 September 2007, P R Truscott, QinetiQ Ltd & Richard Holmberg
//                    Updated extensively prior to this date to deal with
//                    concaved tessellated surfaces, based on the algorithm
//                    of Richard Holmberg.  This had been slightly modified
//                    to determine with inside the geometry by projecting
//                    random rays from the point provided.  Now random rays
//                    are predefined rather than making use of random
//                    number generator at run-time.
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

  SetRandomVectorSet();
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

  SetRandomVectorSet();
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
    solidClosed(false), dirTolerance(1.0E-14)
{
  SetRandomVectorSet();
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
  : G4VSolid(s), fpPolyhedron(0)
{
  dirTolerance = 1.0E-14;
  
  geometryType = "G4TessellatedSolid";
  facets.clear();
  solidClosed  = false;

  cubicVolume = s.cubicVolume;  
  surfaceArea = s.surfaceArea;  

  xMinExtent =  kInfinity;
  xMaxExtent = -kInfinity;
  yMinExtent =  kInfinity;
  yMaxExtent = -kInfinity;
  zMinExtent =  kInfinity;
  zMaxExtent = -kInfinity;

  SetRandomVectorSet();

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
  
  // Copy base class data
  //
  G4VSolid::operator=(s);

  // Copy data
  //
  cubicVolume = s.cubicVolume;  
  surfaceArea = s.surfaceArea;
  fpPolyhedron = 0; 

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
  for (size_t i=0; i<n; i++)
  {
    G4VFacet *facetClone = (s.GetFacet(i))->GetClone();
    AddFacet(facetClone);
  }
  
  if ( s.GetSolidClosed() )  { SetSolidClosed(true); }
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
    aFacet->StreamInfo(G4cout);
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
//
// First the simple test - check if we're outside of the X-Y-Z extremes
// of the tessellated solid.
//
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
#if G4SPECSDEBUG
  G4int nTry                = 7;
#else
  G4int nTry                = 3;
#endif
  G4double distOut          = kInfinity;
  G4double distIn           = kInfinity;
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
  G4int m                   = 0;

  for (G4int i=0; i<nTry; i++)
  {
    G4bool nearParallel = false;
    do
    {
//
//
// We loop until we find direction where the vector is not nearly parallel
// to the surface of any facet since this causes ambiguities.  The usual
// case is that the angles should be sufficiently different, but there are 20
// random directions to select from - hopefully sufficient.
//
      distOut          = kInfinity;
      distIn           = kInfinity;
      G4ThreeVector v  = randir[m];
      m++;
      FacetCI f = facets.begin();
      do
      {
//
//
// Here we loop through the facets to find out if there is an intersection
// between the ray and that facet.  The test if performed separately whether
// the ray is entering the facet or exiting.
//
        crossingO =  ((*f)->Intersect(p,v,true,distO,distFromSurfaceO,normalO));
        crossingI =  ((*f)->Intersect(p,v,false,distI,distFromSurfaceI,normalI));
        if (crossingO || crossingI)
        {
          nearParallel = (crossingO && std::fabs(normalO.dot(v))<dirTolerance) ||
                         (crossingI && std::fabs(normalI.dot(v))<dirTolerance);
          if (!nearParallel)
          {
            if (crossingO && distO > 0.0 && distO < distOut) distOut = distO;
            if (crossingI && distI > 0.0 && distI < distIn)  distIn  = distI;
          }
        }
      } while (!nearParallel && ++f!=facets.end());
    } while (nearParallel && m!=maxTries);

#ifdef G4VERBOSE
    if (m == maxTries)
    {
//
//
// We've run out of random vector directions.  If nTries is set sufficiently
// low (nTries <= 0.5*maxTries) then this would indicate that there is
// something wrong with geometry.
//
      std::ostringstream message;
      G4int oldprc = message.precision(16);
      message << "Cannot determine whether point is inside or outside volume!"
              << G4endl
              << "Solid name       = " << GetName()  << G4endl
              << "Geometry Type    = " << geometryType  << G4endl
              << "Number of facets = " << facets.size() << G4endl
              << "Position:"  << G4endl << G4endl
              << "p.x() = "   << p.x()/mm << " mm" << G4endl
              << "p.y() = "   << p.y()/mm << " mm" << G4endl
              << "p.z() = "   << p.z()/mm << " mm";
      message.precision(oldprc);
      G4Exception("G4TessellatedSolid::Inside()",
                  "InvalidCondition", JustWarning, G4String(message.str()));
    }
#endif
//
//
// In the next if-then-elseif string the logic is as follows:
// (1) You don't hit anything so cannot be inside volume, provided volume
//     constructed correctly!
// (2) Distance to inside (ie. nearest facet such that you enter facet) is
//     shorter than distance to outside (nearest facet such that you exit
//     facet) - on condition of safety distance - therefore we're outside.
// (3) Distance to outside is shorter than distance to inside therefore we're
//     inside.
//
    if (distIn == kInfinity && distOut == kInfinity)
      locationprime = kOutside;
    else if (distIn <= distOut - kCarTolerance*0.5)
      locationprime = kOutside;
    else if (distOut <= distIn - kCarTolerance*0.5)
      locationprime = kInside;

    if (i == 0)  { location = locationprime; }
#ifdef G4VERBOSE
    else if (locationprime != location)
    {
//
//
// Different ray directions result in different answer.  Seems like the
// geometry is not constructed correctly.
//
      std::ostringstream message;
      G4int oldprc = message.precision(16);
      message << "Cannot determine whether point is inside or outside volume!"
              << G4endl
              << "Solid name       = " << GetName()  << G4endl
              << "Geometry Type    = " << geometryType  << G4endl
              << "Number of facets = " << facets.size() << G4endl
              << "Position:"  << G4endl << G4endl
              << "p.x() = "   << p.x()/mm << " mm" << G4endl
              << "p.y() = "   << p.y()/mm << " mm" << G4endl
              << "p.z() = "   << p.z()/mm << " mm";
      message.precision(oldprc);
      G4Exception("G4TessellatedSolid::Inside()",
                  "InvalidCondition", JustWarning, G4String(message.str()));
    }
#endif
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
    std::ostringstream message;
    message << "Point p is not on surface !?" << G4endl
            << "          No facets found for point: " << p << " !" << G4endl
            << "          Returning approximated value for normal.";
    G4Exception("G4TessellatedSolid::SurfaceNormal(p)", "InvalidSetup",
                JustWarning, G4String(message.str()));
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
     std::ostringstream message;
     G4int oldprc = message.precision(16) ;
     message << "Point p is already inside!?" << G4endl
             << "Position:"  << G4endl << G4endl
             << "   p.x() = "   << p.x()/mm << " mm" << G4endl
             << "   p.y() = "   << p.y()/mm << " mm" << G4endl
             << "   p.z() = "   << p.z()/mm << " mm" << G4endl
             << "DistanceToOut(p) == " << DistanceToOut(p);
     message.precision(oldprc) ;
     G4Exception("G4TriangularFacet::DistanceToIn(p,v)", "InvalidSetup",
                 JustWarning, G4String(message.str()));
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    if ((*f)->Intersect(p,v,false,dist,distFromSurface,normal))
    {
//
//
// Set minDist to the new distance to current facet if distFromSurface is in
// positive direction and point is not at surface.  If the point is within
// 0.5*kCarTolerance of the surface, then force distance to be zero and
// leave member function immediately (for efficiency), as proposed by & credit
// to Akira Okumura.
//
      if (distFromSurface > 0.5*kCarTolerance && dist >= 0.0 && dist < minDist)
      {
        minDist  = dist;
      }
      else if (-0.5*kCarTolerance <= dist && dist <= 0.5*kCarTolerance)
      {
        return 0.0;
      }
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
     std::ostringstream message;
     G4int oldprc = message.precision(16) ;
     message << "Point p is already inside!?" << G4endl
             << "Position:"  << G4endl << G4endl
             << "p.x() = "   << p.x()/mm << " mm" << G4endl
             << "p.y() = "   << p.y()/mm << " mm" << G4endl
             << "p.z() = "   << p.z()/mm << " mm" << G4endl
             << "DistanceToOut(p) == " << DistanceToOut(p);
     message.precision(oldprc) ;
     G4Exception("G4TriangularFacet::DistanceToIn(p)", "InvalidSetup",
                 JustWarning, G4String(message.str()));
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    //dist = (*f)->Distance(p,minDist,false);
    dist = (*f)->Distance(p,minDist);
    if (dist < minDist)  { minDist  = dist; }
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
     std::ostringstream message;
     G4int oldprc = message.precision(16) ;
     message << "Point p is already outside!?" << G4endl
             << "Position:"  << G4endl << G4endl
             << "   p.x() = "   << p.x()/mm << " mm" << G4endl
             << "   p.y() = "   << p.y()/mm << " mm" << G4endl
             << "   p.z() = "   << p.z()/mm << " mm" << G4endl
             << "DistanceToIn(p) == " << DistanceToIn(p);
     message.precision(oldprc) ;
     G4Exception("G4TriangularFacet::DistanceToOut(p)", "InvalidSetup",
                 JustWarning, G4String(message.str()));
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
        if (calcNorm) {
          *validNorm = extremeFacets.count(*f);
          *n         = SurfaceNormal(p);
        }  
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
     std::ostringstream message;
     G4int oldprc = message.precision(16) ;
     message << "Point p is already outside!?" << G4endl
             << "Position:"  << G4endl << G4endl
             << "p.x() = "   << p.x()/mm << " mm" << G4endl
             << "p.y() = "   << p.y()/mm << " mm" << G4endl
             << "p.z() = "   << p.z()/mm << " mm" << G4endl
             << "DistanceToIn(p) == " << DistanceToIn(p);
     message.precision(oldprc) ;
     G4Exception("G4TriangularFacet::DistanceToOut(p)", "InvalidSetup",
                 JustWarning, G4String(message.str()));
  }
#endif

  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    //dist = (*f)->Distance(p,minDist,true);
    dist = (*f)->Distance(p,minDist);
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
// Make a clone of the object
//
G4VSolid* G4TessellatedSolid::Clone() const
{
  return new G4TessellatedSolid(*this);
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
  polyhedron->SetReferences();  
 
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
// CalculateExtent
//
// Based on correction provided by Stan Seibert, University of Texas.
//
G4bool
G4TessellatedSolid::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimit,
                                    const G4AffineTransform& pTransform,
                                          G4double& pMin, G4double& pMax) const
{
    G4ThreeVectorList transVertexList(vertexList);

    // Put solid into transformed frame
    for (size_t i=0; i<vertexList.size(); i++)
      { pTransform.ApplyPointTransform(transVertexList[i]); }

    // Find min and max extent in each dimension
    G4ThreeVector minExtent(kInfinity, kInfinity, kInfinity);
    G4ThreeVector maxExtent(-kInfinity, -kInfinity, -kInfinity);
    for (size_t i=0; i<transVertexList.size(); i++)
    {
      for (G4int axis=G4ThreeVector::X; axis < G4ThreeVector::SIZE; axis++)
      {
        G4double coordinate = transVertexList[i][axis];
        if (coordinate < minExtent[axis])
          { minExtent[axis] = coordinate; }
        if (coordinate > maxExtent[axis])
          { maxExtent[axis] = coordinate; }
      }
    }

    // Check for containment and clamp to voxel boundaries
    for (G4int axis=G4ThreeVector::X; axis < G4ThreeVector::SIZE; axis++)
    {
      EAxis geomAxis = kXAxis; // G4 geom classes use different index type
      switch(axis)
      {
        case G4ThreeVector::X: geomAxis = kXAxis; break;
        case G4ThreeVector::Y: geomAxis = kYAxis; break;
        case G4ThreeVector::Z: geomAxis = kZAxis; break;
      }
      G4bool isLimited = pVoxelLimit.IsLimited(geomAxis);
      G4double voxelMinExtent = pVoxelLimit.GetMinExtent(geomAxis);
      G4double voxelMaxExtent = pVoxelLimit.GetMaxExtent(geomAxis);

      if (isLimited)
      {
        if ( minExtent[axis] > voxelMaxExtent+kCarTolerance ||
             maxExtent[axis] < voxelMinExtent-kCarTolerance    )
        {
          return false ;
        }
        else
        {
          if (minExtent[axis] < voxelMinExtent)
          {
            minExtent[axis] = voxelMinExtent ;
          }
          if (maxExtent[axis] > voxelMaxExtent)
          {
            maxExtent[axis] = voxelMaxExtent;
          }
        }
      }
    }

    // Convert pAxis into G4ThreeVector index
    G4int vecAxis=0;
    switch(pAxis)
    {
      case kXAxis: vecAxis = G4ThreeVector::X; break;
      case kYAxis: vecAxis = G4ThreeVector::Y; break;
      case kZAxis: vecAxis = G4ThreeVector::Z; break;
      default: break;
    } 

    pMin = minExtent[vecAxis] - kCarTolerance;
    pMax = maxExtent[vecAxis] + kCarTolerance;

    return true;
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

  G4int i = G4RandFlat::shootInt(facets.size());
  return facets[i]->GetPointOnFace();
}
///////////////////////////////////////////////////////////////////////////////
//
// SetRandomVectorSet
//
// This is a set of predefined random vectors (if that isn't a contradition
// in terms!) used to generate rays from a user-defined point.  The member
// function Inside uses these to determine whether the point is inside or
// outside of the tessellated solid.  All vectors should be unit vectors.
//
void G4TessellatedSolid::SetRandomVectorSet()
{
  randir[0]  = G4ThreeVector(-0.9577428892113370, 0.2732676269591740, 0.0897405271949221);
  randir[1]  = G4ThreeVector(-0.8331264504940770,-0.5162067214954600,-0.1985722492445700);
  randir[2]  = G4ThreeVector(-0.1516671651108820, 0.9666292616127460, 0.2064580868390110);
  randir[3]  = G4ThreeVector( 0.6570250350323190,-0.6944539025883300, 0.2933460081893360);
  randir[4]  = G4ThreeVector(-0.4820456281280320,-0.6331060000098690,-0.6056474264406270);
  randir[5]  = G4ThreeVector( 0.7629032554236800, 0.1016854697539910,-0.6384658864065180);
  randir[6]  = G4ThreeVector( 0.7689540409061150, 0.5034929891988220, 0.3939600142169160);
  randir[7]  = G4ThreeVector( 0.5765188359255740, 0.5997271636278330,-0.5549354566343150);
  randir[8]  = G4ThreeVector( 0.6660632777862070,-0.6362809868288380, 0.3892379937580790);
  randir[9]  = G4ThreeVector( 0.3824415020414780, 0.6541792713761380,-0.6525243125110690);
  randir[10] = G4ThreeVector(-0.5107726564526760, 0.6020905056811610, 0.6136760679616570);
  randir[11] = G4ThreeVector( 0.7459135439578050, 0.6618796061649330, 0.0743530220183488);
  randir[12] = G4ThreeVector( 0.1536405855311580, 0.8117477913978260,-0.5634359711967240);
  randir[13] = G4ThreeVector( 0.0744395301705579,-0.8707110101772920,-0.4861286795736560);
  randir[14] = G4ThreeVector(-0.1665874645185400, 0.6018553940549240,-0.7810369397872780);
  randir[15] = G4ThreeVector( 0.7766902003633100, 0.6014617505959970,-0.1870724331097450);
  randir[16] = G4ThreeVector(-0.8710128685847430,-0.1434320216603030,-0.4698551243971010);
  randir[17] = G4ThreeVector( 0.8901082092766820,-0.4388411398893870, 0.1229871120030100);
  randir[18] = G4ThreeVector(-0.6430417431544370,-0.3295938228697690, 0.6912779675984150);
  randir[19] = G4ThreeVector( 0.6331124368380410, 0.6306211461665000, 0.4488714875425340);

  maxTries = 20;
}
