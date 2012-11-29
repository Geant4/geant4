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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4SurfaceBoundary.cc
//
// ----------------------------------------------------------------------

#include "G4SurfaceBoundary.hh"
#include "geomdefs.hh"
#include "G4CompositeCurve.hh"


G4SurfaceBoundary::G4SurfaceBoundary()
{
}

G4SurfaceBoundary::~G4SurfaceBoundary()
{
}

void G4SurfaceBoundary::Init(const G4CurveVector& bounds0)
{
  bounds= bounds0;
  lastIntersection.Reset();
  
  const G4BoundingBox3D* b= bounds[0]->BBox();
  bBox.Init(b->GetBoxMin(), b->GetBoxMax());
  
  size_t i;
  for ( i=1; i<bounds.size(); i++) 
  {
    b= bounds[i]->BBox();
    bBox.Extend(b->GetBoxMin());
    bBox.Extend(b->GetBoxMax());
  }

  // the points array is probably unused, so the following code is useless
  G4int cnt= 0;

  size_t entr = bounds.size();

  for (i=0; i < entr; i++) 
  {
    G4Curve* c = bounds[i];

    if (c->GetEntityType() == "G4CompositeCurve") 
    {
      G4CompositeCurve* cc = (G4CompositeCurve*)c;
      const G4CurveVector& segments = cc->GetSegments();
      cnt+= segments.size();
    } 
    else 
      cnt++;
  }

  points.resize(cnt);
  
  G4int j= 0;
  
  for (i=0; i<bounds.size(); i++) 
  {
    G4Curve* c= bounds[i];
    if (c->GetEntityType() == "G4CompositeCurve") 
    {
      G4CompositeCurve* cc = (G4CompositeCurve*)c;
      const G4CurveVector& segments = cc->GetSegments();
     
      for (size_t k=0; k<segments.size(); k++) 
      {
	G4Curve* ccc = segments[k];
	G4Point3D p  = ccc->GetEnd();
	points[j]= p;
	j++;
      }

    } 
    else
    {
      G4Point3D p= c->GetEnd();
      points[j]= p;
      j++;
    }
  }
}


G4SurfaceBoundary* G4SurfaceBoundary::Project(const G4Transform3D& tr)
{
  G4CurveVector newBounds;
  G4Curve* a = 0;
  G4Curve* c = 0;
  
  for (size_t i=0; i<bounds.size(); i++)
  {
    c= bounds[i]->Project(tr);
    
    if (c==0) 
    {
      // Remove newBounds and delete all its contents
      while (newBounds.size()>0)
      {
        a = newBounds.back();
        newBounds.pop_back();
        for (G4CurveVector::iterator it=newBounds.begin();
                                     it!=newBounds.end();)
        {
          if (*it==a)
          {
	    it = newBounds.erase(it);
          }
          else
          {
	    it--;
          }
        } 
        if ( a )  { delete a; }
      } 
      return 0;
    }
    // L. Broglia
    c->SetSameSense( bounds[i]->GetSameSense() );

    newBounds.push_back(c);
  }
  
  G4SurfaceBoundary* lof= new G4SurfaceBoundary;
  lof->Init(newBounds);
  return lof;
}

/*
void G4SurfaceBoundary::IntersectRay2D(const G4Ray& ray,
				       G4CurveRayIntersection& is)
{
  is.Reset();
  G4int entr = bounds.entries();
  for (G4int i=0; i < entr; i++) 
  {
    G4Curve& c = *bounds.at(i);
    G4CurveRayIntersection isTmp(c, ray);
    c.IntersectRay2D(ray, isTmp);
    
    if (std::fabs(isTmp.GetDistance()) < std::fabs(is.GetDistance())) 
      is= isTmp;
  }

  lastIntersection= is;
}
*/

G4int G4SurfaceBoundary::IntersectRay2D(const G4Ray& ray)
{
  G4int nbinter = 0;
  G4int temp = 0;

  for (size_t i=0; i < bounds.size(); i++) 
  {   
    G4Curve& c = *bounds[i];
    temp = c.IntersectRay2D(ray);

    // test if the point is on the boundary
    if ( temp == 999 )
      nbinter = 1;
    else
      nbinter += temp; 
  }

  return nbinter;
}


G4bool G4SurfaceBoundary::Tangent(G4CurvePoint&, G4Vector3D& v)
{
  if (lastIntersection.GetDistance() == kInfinity) 
    return false;
  
  return lastIntersection.GetCurve().Tangent(lastIntersection, v);
  // should be true
  
  // cp is ignored for the moment
}


void G4SurfaceBoundary::SplitWithPlane(const G4Point3D&,
				       const G4Vector3D&,
				       G4SurfaceBoundary*&,
				       G4SurfaceBoundary*&)
{
  G4Exception("G4SurfaceBoundary::SplitWithPlane()", "GeomSolids0001",
              FatalException, "Sorry, not yet implemented.");
}

void G4SurfaceBoundary::SplitWithCylinder(const G4CylindricalSurface&,
                                          G4SurfaceBoundary*&, 
                                          G4SurfaceBoundary*&)
{
  G4Exception("G4SurfaceBoundary::SplitWithCylinder()", "GeomSolids0001",
              FatalException, "Sorry, not yet implemented.");
}
