#include "G4SurfaceBoundary.hh"
#include "geomdefs.hh"
#include "G4CompositeCurve.hh"


G4SurfaceBoundary::G4SurfaceBoundary(){}

void G4SurfaceBoundary::Init(const G4CurveVector& bounds0)
{
  bounds= bounds0;
  lastIntersection.Reset();
  
  const G4BoundingBox3D* b= bounds[0]->BBox();
  bBox.Init(b->GetBoxMin(), b->GetBoxMax());
  
  G4int i;
  for ( i=1; i<bounds.entries(); i++) 
  {
    b= bounds[i]->BBox();
    bBox.Extend(b->GetBoxMin());
    bBox.Extend(b->GetBoxMax());
  }

  // the points array is probably unused, so the following code is useless
  G4int cnt= 0;

  G4int entr = bounds.entries();

  for (i=0; i < entr; i++) 
  {
    G4Curve* c = bounds(i);

    if (c->GetEntityType() == "G4CompositeCurve") 
    {
      G4CompositeCurve* cc = (G4CompositeCurve*)c;
      const G4CurveVector& segments = cc->GetSegments();
      cnt+= segments.entries();
    } 
    else 
      cnt++;
  }

  points.reshape(cnt);
  
  G4int j= 0;
  
  for (i=0; i<bounds.entries(); i++) 
  {
    G4Curve* c= bounds(i);
    if (c->GetEntityType() == "G4CompositeCurve") 
    {
      G4CompositeCurve* cc = (G4CompositeCurve*)c;
      const G4CurveVector& segments = cc->GetSegments();
     
      for (G4int i=0; i<segments.entries(); i++) 
      {
	G4Curve* ccc = segments(i);
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


G4SurfaceBoundary::~G4SurfaceBoundary() {}


G4SurfaceBoundary* G4SurfaceBoundary::Project(const G4Transform3D& tr)
{
  G4CurveVector newBounds;
  
  for (G4int i=0; i<bounds.entries(); i++)
  {
    G4Curve* c= bounds[i]->Project(tr);
    
    if (c==0) 
    {
      newBounds.clearAndDestroy();
      return 0;
    }
    // L. Broglia
    c->SetSameSense( bounds[i]->GetSameSense() );

    newBounds.insert(c);
  }
  
  G4SurfaceBoundary* lof= new G4SurfaceBoundary;
  lof->Init(newBounds);
  return lof;
}


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
    
    if (fabs(isTmp.GetDistance()) < fabs(is.GetDistance())) 
      is= isTmp;
  }

  lastIntersection= is;
}


G4bool G4SurfaceBoundary::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  if (lastIntersection.GetDistance() == kInfinity) 
    return false;
  
  return lastIntersection.GetCurve().Tangent(lastIntersection, v);
  // should be true
  
  // cp is ignored for the moment
}


void G4SurfaceBoundary::SplitWithPlane(const G4Point3D& p0,
				       const G4Vector3D& n,
				       G4SurfaceBoundary*& new1,
				       G4SurfaceBoundary*& new2)
{
  G4Exception("To be implemented");
}


