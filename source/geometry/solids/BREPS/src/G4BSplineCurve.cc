// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BSplineCurve.cc,v 1.5 2000-11-08 14:22:08 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BSplineCurve.cc
//
// ----------------------------------------------------------------------

#include "G4BSplineCurve.hh"
#include "G4ControlPoints.hh"
#include "G4KnotVector.hh"

G4BSplineCurve::G4BSplineCurve()
  : degree(0), controlPointsList(0), knots(0), weightsData(0)
{
}

void G4BSplineCurve::Init(G4int degree0, G4Point3DVector* controlPointsList0,
			  G4doubleVector* knots0,
			  G4doubleVector* weightsData0)
{
  degree= degree0;
   
  G4int nbpoints =  controlPointsList0->length();
  controlPointsList = new G4Point3DVector(nbpoints);

  G4int a;  
  for(a = 0; a < nbpoints; a++)
  {
    (*controlPointsList)[a] = (*controlPointsList0)[a];
  }
 
  G4int nbknots = knots0->length();
  knots = new G4doubleVector(nbknots);
  for(a = 0; a < nbknots; a++)
  {
    (*knots)[a] = (*knots0)[a];
  }

  G4int nbweights = weightsData0->length();
  weightsData  = new G4doubleVector(nbweights);
  for(a = 0; a < nbweights; a++)
  {
    (*weightsData)[a] = (*weightsData0)[a];
  }
  

  SetBounds((*knots)[0], (*knots)[knots->length()-1]);
}


G4BSplineCurve::~G4BSplineCurve()
{
  delete [] controlPointsList;
  delete [] knots;
  delete [] weightsData;
}


G4BSplineCurve::G4BSplineCurve(const G4BSplineCurve& right)
{
  delete [] controlPointsList;
  delete [] knots;
  delete [] weightsData;
  Init(right.degree, right.controlPointsList,
       right.knots, right.weightsData);
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;
}

G4BSplineCurve& G4BSplineCurve::operator=(const G4BSplineCurve& right)
{
  if (&right == this) return *this;
  delete [] controlPointsList;
  delete [] knots;
  delete [] weightsData;
  Init(right.degree, right.controlPointsList,
       right.knots, right.weightsData);
  bBox      = right.bBox;
  start     = right.start;
  end       = right.end;
  pStart    = right.pStart;
  pEnd      = right.pEnd;
  pRange    = right.pRange;
  bounded   = right.bounded;
  sameSense = right.sameSense;

  return *this;
}

// add by L. Broglia to pass linkage

G4double G4BSplineCurve::GetPMax() const
{
  return 0.0;
}

G4Point3D G4BSplineCurve::GetPoint(G4double param) const
{
  return G4Point3D(0, 0, 0);
}

G4double  G4BSplineCurve::GetPPoint(const G4Point3D& p) const
{
  return 0.0;
}

/*
#include "G4CurveRayIntersection.hh"

void G4BSplineCurve::IntersectRay2D(const G4Ray& ray,
                                    G4CurveRayIntersection& is)
{
}
*/

G4int G4BSplineCurve::IntersectRay2D(const G4Ray& ray)
{
  // L. Broglia
  G4cout<<"\nWarning ! G4BSplineCurve::IntersectRay2D is empty.";
  return 0;
}

/*
void G4BSplineCurve::CalcCurvePlaneNormal()
{
        //Calc Normal for surface which is used for the projection
    G4ThreeVec norm;
    G4Point3d Pt1 = ControlPointList->get(0,0);
    G4Point3d Pt2 = ControlPointList->get(0,1);
    G4Point3d Pt3 = ControlPointList->get(0,2);    
    G4Point3d a(Pt2.X()-Pt1.X(), Pt2.Y()-Pt1.Y(), Pt2.Z()-Pt1.Z());
    G4Point3d b(Pt3.X()-Pt1.X(), Pt3.Y()-Pt1.Y(), Pt3.Z()-Pt1.Z()); 
    norm.X((a.Y()*b.Z() - a.Z()*b.Y()));
    norm.Y((a.X()*b.Z() - a.Z()*b.X()));
    norm.Z((a.X()*b.Y() - a.Y()*b.X()));

}
*/


G4Curve* G4BSplineCurve::Project(const G4Transform3D& tr)
{
  // just transform + project all control points
  // what about self intersections?
  
  G4int            n                    = controlPointsList->length();
  G4Point3DVector* newControlPointsList = new G4Point3DVector(n);

  for (G4int i=0; i<n; i++) 
  {
    G4Point3D& p= (*newControlPointsList)(i);
    p= tr*(*controlPointsList)(i);
    p.setZ(0);
  }

  G4doubleVector* newKnots= new G4doubleVector(*knots);
  G4doubleVector* newWeightsData= 
    weightsData ? new G4doubleVector(*weightsData) : 0;

  G4BSplineCurve* r= new G4BSplineCurve;
  r->Init(degree, newControlPointsList, newKnots, newWeightsData);

  if (IsBounded()) 
  {
    r->SetBounds(GetPStart(), GetPEnd());
  }
  return r;
}


/*
void G4BSplineCurve::ProjectCurve(const G4Plane& Pl1, const G4Plane& Pl2)
{
    int rows  = ControlPointList->GetRows();
    int cols  = ControlPointList->GetCols();
    int NumberOfPoints = cols * rows;
    ProjectedControlPoints = new G4Point2d*[NumberOfPoints];
    // Loop through points and do projection
    for(int a = 0; a<NumberOfPoints;a++)
    {
	// Create 2d-point
	ProjectedControlPoints[a] = new G4Point2d;
	// Project 3d points into 2d
	Project((*ProjectedControlPoints[a]), ControlPointList->get(0,a), Pl1, Pl2); 
    }
}


int G4BSplineCurve::Inside( const G4Point3d& Hit, const G4Ray& rayref)
{
  const G4Plane& Pl1 = rayref.GetPlane(0);
  const G4Plane& Pl2 = rayref.GetPlane(1);    
  register G4double DistA1, DistA2, DistB1, DistB2;
  // Calc distance from Start point to ray planes
  DistA1 = Start.PlaneDistance(Pl1);
  // Calc distance from End point to ray planes    
  DistB1 = End.PlaneDistance(Pl1);
  if((DistA1<0 && DistB1>0)||(DistA1>0 && DistB1 <0))
    {
      DistA2 = Start.PlaneDistance(Pl2);    
      DistB2 = End.PlaneDistance(Pl2);          
      // This checks the line Start-End of the convex hull
	if(DistA2<0&&DistB2<0)
	  return 1;
    }
  // Test for the other lines of the convex hull
  // If one of them is on a different side than the
  // previously checked line, the curve has to be evaluated
  // against the G4Plane.
  int Points = ControlPointList->GetCols();

  G4Point *CPoint1, *CPoint2;    

  register G4double CDistA1,CDistA2, CDistB1, CDistB2;
  int Flag=0;
  for(int a=0;a<Points-1;a++) 
    {
      CPoint1 = &ControlPointList->get(0,a);
      CPoint2 = &ControlPointList->get(0,a+1);
      CDistA1 = CPoint1->PlaneDistance(Pl1);
      CDistB1 = CPoint2->PlaneDistance(Pl1);
      if((CDistA1<0 && CDistB1>0)||(CDistA1>0 && CDistB1<0))
	{
	  CDistA2 = CPoint1->PlaneDistance(Pl2);
	  CDistB2 = CPoint2->PlaneDistance(Pl2);      
	  if (!(CDistA2<0&&CDistB2<0))
	  {
	    Flag=1;
	    break;
	  }
	}
    }
    	

  if(!Flag)
    return 1;  
  else
    {
      // Evaluate curve & Pl1 intersection, Calc the intersections distance
      // from Pl2 to check which side it lies on.
	G4Point3d IntPoint;
      //      G4cout << "\nG4B_SplineCurve.cc:Inside - Evaluation not yet implemented!!!\n";
      // IntPoint = ...
      G4double IntDist =  IntPoint.PlaneDistance(Pl2);
      if(IntDist<0)
	return 1;
    }
  return 0;
    
}
*/


void G4BSplineCurve::InitBounded()
{
  // just like in the old functions
  G4int pointCount = controlPointsList->length();
  bBox.Init( (*controlPointsList)(0) );
  for (G4int i=1; i<pointCount; i++) 
  {
    bBox.Extend( (*controlPointsList)(i) );
  }
}


/*
G4Point3d G4BSplineCurve::GetBoundMin()
{
  G4Point3d Min = PINFINITY;
  int PointCount = ControlPointList->GetCols();
  G4Point3d Tmp;
  for(int a=0;a<PointCount;a++)
    {
	Tmp = ControlPointList->get(0,a);
	Min > Tmp;
    }
  return Min;
}

G4Point3d G4BSplineCurve::GetBoundMax()
{
  G4Point3d Max = -PINFINITY;
  G4Point3d Tmp;
  int PointCount = ControlPointList->GetCols();
  for(int a=0;a<PointCount;a++)
    {
	Tmp = ControlPointList->get(0,a);
	Max > Tmp;
    }
  return Max;
}
*/


G4bool G4BSplineCurve::Tangent(G4CurvePoint& cp, G4Vector3D& v)
{
  G4Exception("G4BSplineCurve::Tangent");
  return false;
}


