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
// * subject DEFCON 705 IPR conditions.                               *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.cc,v 1.4 2006/06/29 18:48:57 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//
// Standard contructor has blank name and defines no facets.
//
G4TessellatedSolid::G4TessellatedSolid ()
  : G4VSolid("dummy"), fpPolyhedron(0), cubicVolume(0.)
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
  : G4VSolid(name), fpPolyhedron(0), cubicVolume(0.)
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
// Destructor.
//
G4TessellatedSolid::~G4TessellatedSolid ()
{
  DeleteObjects ();
}

///////////////////////////////////////////////////////////////////////////////
//
// Define assignment operator.
//
const G4TessellatedSolid &G4TessellatedSolid::operator=
   (const G4TessellatedSolid &s)
{
  if (&s == this) return *this;
  
  DeleteObjects ();
  CopyObjects (s);
  
  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
void G4TessellatedSolid::DeleteObjects ()
{
  for (std::vector<G4VFacet *>::iterator f=facets.begin(); 
       f!=facets.end(); f++) delete *f;
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
    solidClosed = true;
  }
  else
  {
    solidClosed = false;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4bool G4TessellatedSolid::GetSolidClosed () const
  {return solidClosed;}

///////////////////////////////////////////////////////////////////////////////
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
G4VFacet *G4TessellatedSolid::GetFacet (size_t i) const
{
  return facets[i];
}

///////////////////////////////////////////////////////////////////////////////
//
size_t G4TessellatedSolid::GetNumberOfFacets () const
{
  return facets.size();
}

///////////////////////////////////////////////////////////////////////////////
//
EInside G4TessellatedSolid::Inside (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double dist    = 0.0;
  typedef std::multimap< G4double, FacetCI, std::less<G4double> > DistMapType;
  DistMapType distmap;
  size_t purgeIntv = 25;
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist);
    distmap.insert(DistMapType::value_type(dist,f));
    minDist = distmap.begin()->first;
    if (distmap.size() > purgeIntv)
    {
      DistMapType::iterator it =
        distmap.lower_bound(minDist + 0.5*kCarTolerance);
      it++;
      if (it != distmap.end())
      {
        DistMapType::iterator itend = distmap.end();
        itend--;
        distmap.erase (it,itend);
      }
      if (distmap.size() > purgeIntv) purgeIntv = 2*distmap.size();
    }
  }

  EInside inside = kInside;
  
  if (minDist <= 0.5*kCarTolerance) {inside = kSurface;}
  else
  {
    DistMapType::const_iterator itcut = 
      distmap.lower_bound(minDist + 0.5* kCarTolerance);
    itcut++;
    DistMapType::const_iterator it = distmap.begin();
    do
    {
      if (!((*(it->second))->IsInside(p))) {inside = kOutside;}
    } while (inside == kInside && ++it != itcut);
  }

  return inside;
}

///////////////////////////////////////////////////////////////////////////////
//
G4ThreeVector G4TessellatedSolid::SurfaceNormal (const G4ThreeVector &p) const
{
  FacetCI minFacet;
  G4double minDist   = kInfinity;
  G4double dist      = 0.0;
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist);
    if (dist < minDist)
    {
      minDist  = dist;
      minFacet = f;
    }
  }
  
  return (*minFacet)->GetSurfaceNormal();
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToIn (const G4ThreeVector &p,
  const G4ThreeVector &v) const
{
  G4double minDist         = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal(0.0,0.0,0.0);
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    if ((*f)->Intersect(p,v,false,dist,distFromSurface,normal))
    {
      if (dist < minDist) minDist  = dist;
    }
  }

  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToIn (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double dist    = 0.0;
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist,false);
    if (dist < minDist) minDist  = dist;
  }
  
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToOut (const G4ThreeVector &p,
                    const G4ThreeVector &v, const G4bool calcNorm,
                          G4bool *validNorm, G4ThreeVector *n) const
{
  G4double minDist1        = kInfinity;
  G4double minDist2        = kInfinity;
  G4double dist            = 0.0;
  G4double distFromSurface = 0.0;
  G4ThreeVector normal(0.0,0.0,0.0);
  G4ThreeVector minNormal1(0.0,0.0,0.0);
  G4ThreeVector minNormal2(0.0,0.0,0.0);
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    if ((*f)->Intersect(p,v,true,dist,distFromSurface,normal))
    {
      if (dist < minDist1)
      {
        if (v.dot(normal) > dirTolerance)
        {
          minDist1   = dist;
          minNormal1 = normal;
        }
        else if (dist < minDist2)
        {
          minDist2   = dist;
          minNormal2 = normal;
        }
      }
    }
  }
  
  if (minDist1 < kInfinity)
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n         = minNormal1;
    }
    return minDist1;
  }
  else
  {
    if (calcNorm)
    {
      *validNorm = true;
      *n         = minNormal2;
    }
    return minDist2;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
G4double G4TessellatedSolid::DistanceToOut (const G4ThreeVector &p) const
{
  G4double minDist = kInfinity;
  G4double dist    = 0.0;
  
  for (FacetCI f=facets.begin(); f!=facets.end(); f++)
  {
    dist = (*f)->Distance(p,minDist,true);
    if (dist < minDist) minDist  = dist;
  }
  
  return minDist;
}

///////////////////////////////////////////////////////////////////////////////
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
//G4double G4TessellatedSolid::GetCubicVolume ()
//  {return cubicVolume;}
  
