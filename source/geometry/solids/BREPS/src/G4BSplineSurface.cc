// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BSplineSurface.cc,v 1.7 2000-02-25 15:58:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4BSplineSurface.hh"
#include "G4BezierSurface.hh"
#include "G4ControlPoints.hh"
#include "G4BoundingBox3D.hh"

G4BSplineSurface::G4BSplineSurface()
{
  distance = kInfinity;
  dir=ROW;
  first_hit = Hit = (G4UVHit*)0;
  ctl_points = (G4ControlPoints*)0;
  u_knots = v_knots = tmp_knots = (G4KnotVector*)0;
}


G4BSplineSurface::G4BSplineSurface(const char* nurbfilename, G4Ray& rayref)
{
  distance = kInfinity;    
  first_hit = Hit = (G4UVHit*)0;
  ctl_points = (G4ControlPoints*)0;
  u_knots = v_knots = tmp_knots = (G4KnotVector*)0;
}


G4BSplineSurface::G4BSplineSurface(const  G4BSplineSurface &tmp) 
{
  distance = tmp.distance;
  first_hit = Hit = (G4UVHit*)0;

  //    next=this;
  order[0] = tmp.order[0];
  order[1] = tmp.order[1];
  dir = tmp.dir;

  u_knots    = new G4KnotVector(*tmp.u_knots);
  v_knots    = new G4KnotVector(*tmp.v_knots);
  tmp_knots  = (G4KnotVector*)0;

  ctl_points = new G4ControlPoints(*tmp.ctl_points);
}


G4BSplineSurface::G4BSplineSurface(G4int u, G4int v, G4KnotVector& u_kv, 
				   G4KnotVector& v_kv, G4ControlPoints& cp)
{
  first_hit = Hit = (G4UVHit*)0;

  order[0] = u+1;
  order[1] = v+1;

  u_knots    = new G4KnotVector(u_kv);
  v_knots    = new G4KnotVector(v_kv);
  tmp_knots  = (G4KnotVector*)0;
  
  ctl_points =  new G4ControlPoints(cp);
}


G4BSplineSurface::~G4BSplineSurface()
{
  delete u_knots;
  delete v_knots;
  delete ctl_points;
  G4UVHit* temphit=Hit;
  Hit = first_hit;
  while(Hit!=(G4UVHit*)0)
  {
    Hit=Hit->next;
    delete temphit;
    temphit=Hit;
  }
  // delete temphit;// remove last
  
}


int G4BSplineSurface::Intersect(const G4Ray& rayref)
{
  Intersected = 1;
  FindIntersections(rayref);
  G4BezierSurface *bez_ptr;
  bezier_list.MoveToFirst();
  distance = kInfinity;
  
  while( bezier_list.index != (G4Surface*)0)
  {
    bez_ptr = (G4BezierSurface*)bezier_list.GetSurface();
    
    if(bez_ptr->Active())
      if(distance > bez_ptr->Distance())
      {
	// Put data from closest bezier to b-spline data struct
	closest_hit = bez_ptr->AveragePoint();
	distance = bez_ptr->Distance();
      }
      else
      {
	// Set other beziers as inactive
	bez_ptr->Active(0);
	    
	// Remove beziers that are not closest
	//  bezier_list.RemoveSurface(bez_ptr);
      }

    bezier_list.Step();
  }
  
  bezier_list.MoveToFirst();
    
  if(bezier_list.number_of_elements)
    return 1;
  else
  {
    active=0;
    return 0;
  }
}


G4Point3D G4BSplineSurface::FinalIntersection()
{
  // Compute the real intersection point.
  G4BezierSurface* bez_ptr;
  while ( bezier_list.number_of_elements > 0 && 
	  bezier_list.index != (G4Surface*)0)
  {
    bez_ptr = (G4BezierSurface*)bezier_list.GetSurface();
    int tmp = 0;

    // L. Broglia
    // Modify G4BezierSurface intersection function name 
    // tmp = bez_ptr->Intersect( bezier_list);
    tmp = bez_ptr->BIntersect( bezier_list);

    if(!tmp)
    {
      bezier_list.RemoveSurface(bez_ptr);
      if(bezier_list.index != (G4Surface*)0)
	bezier_list.index->Active(1);
    }
    else
      if(tmp==1)
      {
	active=1;
	// Hit found
	AddHit(bez_ptr->GetU(), bez_ptr->GetV());

	// Delete beziers
	bezier_list.EmptyList();
      }
      else
	if(tmp==2)
	{	
	  // The bezier was split so the last
	  // two surfaces in the List should
	  // be bbox tested and if passed
	  // clipped in both dirs.
		
	  // Move to first
	  bezier_list.MoveToFirst();
	  // Find the second last.
// What!?  Casting a G4Surface* to a G4SurfaceList*  !?!?!? - GC
//
//	  if(bezier_list.index != bezier_list.last)
//	    while ( ((G4SurfaceList*)bezier_list.index)->next !=
//		    bezier_list.last)  bezier_list.Step();
//
// Try the following instead (if that's the wished behavior)...
//
	  if(bezier_list.index != bezier_list.last)
	    while (bezier_list.next != bezier_list.last)
	      bezier_list.Step();
	  
	  G4BezierSurface* tmp = (G4BezierSurface*) bezier_list.GetSurface();
	  tmp->CalcBBox();
	  
// L. Broglia		tmp->bbox->Test();

	  int result=0;
	  if(tmp->bbox->GetTestResult())
	  {
	    // Clip
	    while(!result)
	      result = tmp->ClipBothDirs();  
	  }
	  else
	  {
	    bezier_list.RemoveSurface(tmp);
	  }
	  // Second surface
	  tmp = (G4BezierSurface*) bezier_list.GetLastSurface();
	  tmp->CalcBBox();

// L. Broglia		tmp->bbox->Test();

	  if(tmp->bbox->GetTestResult())
	  {
	    result = 0;
	    while(!result)
	      result = tmp->ClipBothDirs();
	  }
	  else
	  {
	    bezier_list.RemoveSurface(tmp);
	  }
	  
	  bezier_list.RemoveSurface(bez_ptr);
	  bezier_list.MoveToFirst();
	}
    
    bezier_list.Step();
  }//While....
  
  Hit = first_hit;
  G4Point3D result;
  if(Hit == (G4UVHit*)0)
    active = 0;
  else
  {
    while(Hit != (G4UVHit*)0)
    {
      // L. Broglia
      // Modify function name
      // result = Evaluate();
      result = BSEvaluate();

      Hit = Hit->next;
    }

    Hit = first_hit;
  }

  return result;
}	


void G4BSplineSurface::CalcBBox()
{
    
  // Finds the bounds of the b-spline surface iow
  // calculates the bounds for a bounding box
  // to the surface. The bounding box is used
  // for a preliminary check of intersection.
  
  register G4Point3D box_min = PINFINITY;
  register G4Point3D box_max =-PINFINITY;        
  
  // Loop to search the whole control point mesh
  // for the minimum and maximum values for x, y and z.
  
  for(register int a = ctl_points->GetRows()-1; a>=0;a--)
    for(register int b = ctl_points->GetCols()-1; b>=0;b--)
    {
      G4Point3D tmp = ctl_points->Get3D(a,b);
      if((box_min.x()) > (tmp.x())) box_min.setX(tmp.x());
      if((box_min.y()) > (tmp.y())) box_min.setY(tmp.y());
      if((box_min.z()) > (tmp.z())) box_min.setZ(tmp.z());
      if((box_max.x()) < (tmp.x())) box_max.setX(tmp.x());
      if((box_max.y()) < (tmp.y())) box_max.setY(tmp.y());
      if((box_max.z()) < (tmp.z())) box_max.setZ(tmp.z()); 
    }
  bbox = new G4BoundingBox3D( box_min, box_max);
}


G4ProjectedSurface* G4BSplineSurface::CopyToProjectedSurface
(const G4Ray& rayref)
{
  G4ProjectedSurface* proj_srf = new G4ProjectedSurface() ;
  proj_srf->PutOrder(0,GetOrder(0));
  proj_srf->PutOrder(1,GetOrder(1));
  proj_srf->dir = dir;

  proj_srf->u_knots     =  new G4KnotVector(*u_knots);
  proj_srf->v_knots     =  new G4KnotVector(*v_knots);
  proj_srf->ctl_points  = new G4ControlPoints
    (2, ctl_points->GetRows(), ctl_points->GetCols());

  const G4Plane& plane1 = rayref.GetPlane(1);
  const G4Plane& plane2 = rayref.GetPlane(2);
  ProjectNURBSurfaceTo2D(plane1, plane2, proj_srf);

  return proj_srf;
}


void G4BSplineSurface::FindIntersections(const G4Ray& rayref)
{
  // Do the projection to 2D
  G4ProjectedSurface* proj_srf = CopyToProjectedSurface(rayref);

  // Put surface in projected List
  projected_list.AddSurface(proj_srf);
  
  // Loop through List of projected surfaces
  while(projected_list.number_of_elements > 0)
  {
    // Get first in List
    proj_srf = (G4ProjectedSurface*)projected_list.GetSurface();
    
    // Create the bounding box for the projected surface.
    proj_srf->CalcBBox();
    
// L. Broglia	proj_srf->bbox->Test();
	
    // Check bbox test result is ok
    if(proj_srf->bbox->GetTestResult())
      // Convert the projected surface to a bezier. Split if necessary.
      proj_srf->ConvertToBezier(projected_list, bezier_list);

    // Remove projected surface
    projected_list.RemoveSurface(proj_srf);
  }
  
  // Loop through the bezier List
  G4BezierSurface* bez_ptr;
  distance = kInfinity;

  while(bezier_list.index != (G4Surface*)0)
  {
    bez_ptr = (G4BezierSurface*)bezier_list.GetSurface();

    // Add a temporary Hit
    AddHit(bez_ptr->UAverage(), bez_ptr->VAverage());
	
    // Evaluate Hit

    // L. Broglia
    // Modify function name
    // bez_ptr->SetAveragePoint(Evaluate());
    bez_ptr->SetAveragePoint(BSEvaluate());
    
    // Calculate distance to ray origin
    bez_ptr->CalcDistance(rayref.GetStart());

    // Put closest to b_splines distance value
    if(bez_ptr->Distance() < distance) distance = bez_ptr->Distance();  
    
    // Remove the temporary Hit
    if (first_hit == Hit) first_hit = (G4UVHit*)0;
    delete Hit;
    Hit = (G4UVHit*)0;
    
    // Move to next in the List
    bezier_list.Step();
  }
  
  bezier_list.MoveToFirst();
  if(bezier_list.number_of_elements == 0)
  {
    active=0; 
    return;
  }
  
  // Check that approx Hit is in direction of ray
  const G4Point3D&   Pt         = rayref.GetStart();
  const G4Vector3D&  Dir        = rayref.GetDir();
  G4Point3D          TestPoint  = (0.00001*Dir) + Pt;
  G4BezierSurface*   Bsrf       = (G4BezierSurface*)bezier_list.GetSurface(0);

  G4Point3D          AveragePoint = Bsrf->AveragePoint(); 
  G4double           TestDistance = TestPoint.distance2(AveragePoint);

  if(TestDistance > distance)
    // Hit behind ray starting point, no intersection.
    active=0;
}


void G4BSplineSurface::AddHit(G4double u, G4double v)
{
  if(Hit == (G4UVHit*)0)
  {
    first_hit = new G4UVHit(u,v);
    first_hit->next = (G4UVHit*)0;
    Hit = first_hit;
  }
  else
  {
    Hit->next = new G4UVHit(u,v);
    Hit = Hit->next;
    Hit->next=(G4UVHit*)0;
  }
}


void G4BSplineSurface::ProjectNURBSurfaceTo2D
                                (const G4Plane& plane1, const G4Plane& plane2,
				 register G4ProjectedSurface* proj_srf)
{
  // Projects the nurb surface so that the z-axis = ray. 
  
  /* L. Broglia
  G4Point* tmp = (G4Point*)&ctl_points->get(0,0);
  */

  G4PointRat tmp = ctl_points->GetRat(0,0);
  int rational = tmp.GetType();// Get the type of control point
  register G4Point3D psrfcoords;
  register int rows = ctl_points->GetRows();
  register int cols = ctl_points->GetCols();
  
  for (register int i=0; i< rows; i++)
    for(register int j=0; j < cols;j++)
    {
      if ( rational==4 ) // 4 coordinates
      {
	G4PointRat& srfcoords = ctl_points->GetRat(i, j);
		
// L. Broglia	
// Changes for new G4PointRat

	// Calculate the x- and y-coordinates for the new 
	// 2-D surface.  
	psrfcoords.setX((  srfcoords.x()  * plane1.a 
			  +srfcoords.y()  * plane1.b
			  +srfcoords.z()  * plane1.c
			  -srfcoords.w()  * plane1.d));
	psrfcoords.setY((  srfcoords.x()  * plane2.a
			  +srfcoords.y()  * plane2.b
			  +srfcoords.z()  * plane2.c
			  -srfcoords.w()  * plane2.d));
	
	proj_srf->ctl_points->put(i,j,psrfcoords);
      }
      else  // 3 coordinates
      {
	G4Point3D srfcoords = ctl_points->Get3D(i, j);
	
	psrfcoords.setX((  srfcoords.x()  * plane1.a 
			  +srfcoords.y()  * plane1.b
			  +srfcoords.z()  * plane1.c
			                  - plane1.d));
	
	psrfcoords.setY((  srfcoords.x()  * plane2.a
			  +srfcoords.y()  * plane2.b
			  +srfcoords.z()  * plane2.c
			                  - plane2.d));
	
	proj_srf->ctl_points->put(i,j,psrfcoords);		    
      }
    }
} 

/* L. Broglia
   Changes for new G4PointRat 
G4Point& G4BSplineSurface::InternalEvalCrv(int i, G4ControlPoints* crv)*/

G4PointRat& G4BSplineSurface::InternalEvalCrv(int i, G4ControlPoints* crv)
{
  if ( ord <= 1 ) 
    return crv->GetRat(i, k_index);
  
  register int j = k_index;
  
  while ( j > (k_index - ord + 1)) 
  {
    register G4double  k1, k2;
    
    k1 = tmp_knots->GetKnot((j + ord - 1));
    k2 = tmp_knots->GetKnot(j); 
    
    if ((abs(k1 - k2)) > kCarTolerance )
    { 
      /* L. Broglia		
      register G4PointRat* pts1 = &crv->get(i,j-1);
      register G4PointRat* pts2 = &crv->get(i,j  );	    
      if(pts1->GetType()==3)
      {
	crv->CalcValues(k1, param, *(G4Point3D*)pts1, k2, *(G4Point3D*)pts2);
	crv->put(0, j, *(G4Point3D*)pts2);			     
      }
      else
      {
	crv->CalcValues(k1, param, *(G4PointRat*)pts1, k2, *(G4PointRat*)pts2);
	crv->put(0, j, *(G4PointRat*)pts2);	      
      }
      */
      register G4PointRat* pts1 = &crv->GetRat(i,j-1);
      register G4PointRat* pts2 = &crv->GetRat(i,j  );
    } 		

    j--;
  } 	

  ord = ord-1;
  return InternalEvalCrv(0, crv); // Recursion
}


G4Point3D  G4BSplineSurface::BSEvaluate()
{
  register int           	i;
  register int              	row_size = ctl_points->GetRows();
  register G4ControlPoints      *diff_curve;
  register G4ControlPoints*     curves;
  G4Point3D                     result;

  /* L. Broglia  
  G4Point* tmp = (G4Point*)&ctl_points->get(0,0);
  */
  
  G4PointRat* tmp = &ctl_points->GetRat(0,0);

  register int point_type = tmp->GetType();
  diff_curve = new G4ControlPoints(point_type, row_size, 1);
  k_index = u_knots->GetKnotIndex(Hit->u, GetOrder(ROW) );
  
  ord = GetOrder(ROW);
  if(k_index==-1)
  {
    delete diff_curve;
    active = 0;
    return result;
  }

  curves=new G4ControlPoints(*ctl_points);
  tmp_knots = u_knots;
  param = Hit->u;
  
  if(point_type == 4)
  {
    for ( i = 0; i < row_size; i++)
    {
      ord = GetOrder(ROW);
      register G4PointRat rtr_pt = (G4PointRat&) InternalEvalCrv(i, curves);
      diff_curve->put(0,i,rtr_pt);
    }

    k_index = v_knots->GetKnotIndex( Hit->v, GetOrder(COL) );
    if(k_index==-1)
    {
      delete diff_curve;
      delete curves;
      active = 0;
      return result;
    }
	
    ord = GetOrder(COL);
    tmp_knots = v_knots;
    param = Hit->v;
	
    // Evaluate the diff_curve...
    // G4PointRat rat_result = (G4PointRat&) InternalEvalCrv(0, diff_curve);
    G4PointRat rat_result(InternalEvalCrv(0, diff_curve));
    
    // Calc the 3D values.
    // L. Broglia
    // Changes for new G4PointRat
    result.setX(rat_result.x()/rat_result.w());
    result.setY(rat_result.y()/rat_result.w());
    result.setZ(rat_result.z()/rat_result.w());
  }
  else
    if(point_type == 3)
    {
      for ( i = 0; i < row_size; i++)
      {
	ord = GetOrder(ROW);
	// G4Point3D rtr_pt  = (G4Point3D&) InternalEvalCrv(i, curves);
	G4Point3D rtr_pt = (InternalEvalCrv(i, curves)).pt();
	diff_curve->put(0,i,rtr_pt);
      }
	
      k_index = v_knots->GetKnotIndex( Hit->v, GetOrder(COL) );
      if(k_index==-1)
      {
	delete diff_curve;
	delete curves;
	active = 0;
	return result;
      }
      
      ord = GetOrder(COL);
      tmp_knots = v_knots;
      param = Hit->v;

      // Evaluate the diff_curve...
      result = (InternalEvalCrv(0, diff_curve)).pt();
    }
  
  delete diff_curve;
  delete curves;
  closest_hit = result;
  return result;
}


G4Point3D G4BSplineSurface::Evaluation(const G4Ray& rayref)
{
  // Delete old UVhits
  G4UVHit* temphit=Hit;
  while(Hit!=(G4UVHit*)0)
  {
    Hit=Hit->next;
    delete temphit;
    temphit=Hit;
  }
  
  // delete temphit;
  
  // Get the real Hit point
  closest_hit = FinalIntersection();
  
  // The following part (commented out) is old bullshit
  // Chech that Hit is not in a void i.e. InnerBoundary.
  //    for(int a=0; a<NumberOfInnerBoundaries;a++)
  //      if(InnerBoundary[a]->Inside(closest_hit, rayref))
  //	{
  //	  Active(0);
  //	  Distance(kInfinity);
  //	  return closest_hit;
  //	}
  return closest_hit;
}


G4double G4BSplineSurface::ClosestDistanceToPoint(const G4Point3D& Pt)
{
  G4double PointDistance=0;
  PointDistance = ctl_points->ClosestDistanceToPoint(Pt);
  return PointDistance;
}








