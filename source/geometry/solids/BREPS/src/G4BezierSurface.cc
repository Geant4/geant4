// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BezierSurface.cc,v 1.4 2000-11-08 14:22:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BezierSurface.cc
//
// ----------------------------------------------------------------------
// History:
// -------
// - Replaced addition of coordinates by addition of 2 points
//   (L. Broglia, 10/10/98)
// ----------------------------------------------------------------------

#include "G4BezierSurface.hh"
#include "G4ConvexHull.hh"

G4double G4BezierSurface::Tolerance=0;
G4int G4BezierSurface::Clips=0;
G4int G4BezierSurface::Splits=0;


G4BezierSurface::G4BezierSurface()
{
  oslo_m     = (G4OsloMatrix*)0;
  new_knots  = (G4KnotVector*)0;
  old_points = (G4ControlPoints*)0;

  u[0]=0; u[1]=0;
  v[0]=0; v[1]=0;
}

G4BezierSurface::~G4BezierSurface()
{
  delete u_knots;
  delete v_knots;
  delete new_knots;
  delete ctl_points;
  delete old_points;
  
  G4OsloMatrix* temp_oslo = oslo_m;
  
  while(oslo_m != (G4OsloMatrix*)0)
  {
    oslo_m = oslo_m->GetNextNode();
    delete temp_oslo;
    temp_oslo = oslo_m;
  }

  delete oslo_m;
  delete bbox;
}

G4BezierSurface::G4BezierSurface(const G4BezierSurface&)
{
}

G4Vector3D G4BezierSurface::SurfaceNormal(const G4Point3D& Pt) const
{
  return G4Vector3D(0,0,0);
}

G4int G4BezierSurface::ClipBothDirs()
{
  dir = ROW;
  ClipSurface(); 
  
  //   G4cout << "\n CLIP BOTH DIRS  1: " << smin << "  " << smax;

  if(smin > 1.0 || smax < 0.0)
  {
    bezier_list->RemoveSurface(this);
    return 1;
  }
  else
    if((smax - smin) > 0.8)
    {
      SplitNURBSurface();
      return 0;
    }
  
  LocalizeClipValues();
  SetValues();
  
  // Other G4Vector3D clipping and testing.
  dir = COL;
  ClipSurface();
  //    G4cout << "\n CLIP BOTH DIRS  2: " << smin << "  " << smax;
    
  if(smin > 1.0 || smax < 0.0)
  {
    bezier_list->RemoveSurface(this);
    return 1;
  }
  else
    if((smax - smin) > 0.8)
    {
      SplitNURBSurface();
      return 0;
    }

  LocalizeClipValues();    
  SetValues();
  CalcAverage();
  return 1;
}


void G4BezierSurface::CalcBBox()
{
  // Finds the bounds of the 2D-projected nurb iow
  // calculates the bounds for a bounding rectangle
  // to the surface. The bounding rectangle is used
  // for a preliminary check of intersection.
  register G4Point3D box_min = G4Point3D(PINFINITY);
  register G4Point3D box_max = G4Point3D(-PINFINITY);
 
    
  // Loop to search the whole control point mesh
  // for the minimum and maximum values for.X() and y.
  for(register G4int a = ctl_points->GetRows()-1; a>=0;a--)
    for(register G4int b = ctl_points->GetCols()-1; b>=0;b--)
    {
/* L. Broglia
      G4Point2d& tmp = (G4Point2d&)ctl_points->get(a,b);
      if((box_min.X()) > (tmp.X())) box_min.X(tmp.X());
      if((box_max.X()) < (tmp.X())) box_max.X(tmp.X());	
      if((box_min.Y()) > (tmp.Y())) box_min.Y(tmp.Y());	
      if((box_max.Y()) < (tmp.Y())) box_max.Y(tmp.Y());	
*/
      G4Point3D tmp = ctl_points->Get3D(a,b);
      if((box_min.x()) > (tmp.x())) box_min.setX(tmp.x());
      if((box_max.x()) < (tmp.x())) box_max.setX(tmp.x());	
      if((box_min.y()) > (tmp.y())) box_min.setY(tmp.y());	
      if((box_max.y()) < (tmp.y())) box_max.setY(tmp.y());    
    }
	
  bbox = new G4BoundingBox3D(box_min, box_max);
}


void G4BezierSurface::CalcAverage()
{
  // Calculate the average point from the average clip-values.
  average_u = (u_min + u_max)/2.0;
  average_v = (v_min + v_max)/2.0;    
}


void G4BezierSurface::CalcDistance(const G4Point3D& ray_start)
{
  // Calculate the distance between the average point and
  // the ray starting point.
  distance = ((((ray_start.x() - average_pt.x())*
		(ray_start.x() - average_pt.x()))+
	       ((ray_start.y() - average_pt.y())*
		(ray_start.y() - average_pt.y()))+
	       ((ray_start.z() - average_pt.z())*
		(ray_start.z() - average_pt.z()))));
}


void G4BezierSurface::SetValues()
{
  if(dir)
  {
    v_min = smin;
    v_max = smax;
  }
  else
  {
    u_min = smin;
    u_max = smax;
  }
}

	
G4int G4BezierSurface::BIntersect(G4SurfaceList& bez_list)
{
  bezier_list = &bez_list;
  G4int clip_regions = 0; // Used for tolerance/efficiency-testing
  
  do
  {
    // Calc bbox
    CalcBBox();

    // Test bbox
/* L. Broglia
    bbox->Test2dBBox();
*/
    // bbox->Test();

    // Check result
    if(!bbox->GetTestResult())
      return 0;
    
    // The first clipping has already been Done
    // previously so we continue by doing the
    // actual clip.
    
    // Cut out the clipped region of the surface
    GetClippedRegionFromSurface();
    clip_regions++;
    
    // Calculate the knot vectors and control points
    // for the clipped surface
    RefineSurface();    

    // Gets the u- and v-bounds for the clipped surface
    u_min = u_knots->GetKnot(0);	
    u_max = u_knots->GetKnot(u_knots->GetSize() - 1);	
    v_min = v_knots->GetKnot(0);	
    v_max = v_knots->GetKnot(v_knots->GetSize() - 1); 
    
    // Choose the G4Vector3D for the next() clipping so that
    // the larger side will be clipped.  
    if( (u_max - u_min) < (v_max - v_min) )	
      dir = 1;
    else
      dir = 0;

    // Calculate the clip points
    ClipSurface();
    //	    G4cout << "\n       SMINMAX : " << smin << "  " << smax; 
    
    // The ray intersects with the bounding box
    // but not with the surface itself.   
    if( smin > 1.0 || smax < 0.0 )
    {
      //	    G4cout << "\nG4BezierSurface::Intersect : bezier missed!"; 
      //	    bezier_list->RemoveSurface(this);
      return 0;
    }
    
    if( (smax - smin) > 0.8)
    {
      // Multiple intersections
      //	    G4cout << "\nG4BezierSurface::Intersect : Bezier split.";
      SplitNURBSurface();
      // Now the two new surfaces should also be
      // clipped in both G4Vector3Ds i.e the
      // last and the second last surface
      // in the List. This is Done after returning
      // from this function.
      //	    G4cout << "\n\n  BEZ SPLIT in final Calc! \n\n";

	    
      return 2;
    }
    
    // Calculate the smin and smax values on the
    // b_spline.
    LocalizeClipValues();	
    
    // Check if the size of the remaining surface is within the
    // Tolerance .  
  } while ((u_max - u_min > Tolerance) || (v_max - v_min) > Tolerance);    
  
  SetValues();
  //    G4cout << "\nG4BezierSurface::Intersect :Regions were cut " 
  //           << clip_regions << "  Times.\n";
  
  return 1;
}


void G4BezierSurface::ClipSurface()
{
  // This routine is described in Computer Graphics, Volume 24, 
  // Number 4, August 1990 under the title Ray Tracing Trimmed
  // Rational Surface Patches.
  

  //    G4cout << "\nBezier clip.";
  
  register G4int i,j;
  register G4ConvexHull *ch_ptr, *ch_tmp, *ch_first;
  register G4int col_size = ctl_points->GetCols();
  register G4int row_size = ctl_points->GetRows();
  
  // The four cornerpoints of the controlpoint mesh.

/* L. Broglia
  register G4Point2d pt1 = ctl_points->get(0,0);    
  register G4Point2d pt2 = ctl_points->get(0,col_size-1);    
  register G4Point2d pt3 = ctl_points->get(row_size-1,0);    
  register G4Point2d pt4 = ctl_points->get(row_size-1,col_size-1);    
  register G4Point2d v1,v2,v3;
*/
  register G4Point3D pt1 = ctl_points->Get3D(0,0);    
  register G4Point3D pt2 = ctl_points->Get3D(0,col_size-1);    
  register G4Point3D pt3 = ctl_points->Get3D(row_size-1,0);    
  register G4Point3D pt4 = ctl_points->Get3D(row_size-1,col_size-1);    
  register G4Point3D v1,v2,v3;

  if ( dir == ROW)
  {
    // Vectors from cornerpoints
    v1 = (pt1 - pt3);
    //	v1.X() = pt1.X() - pt3.X();
    //	v1.Y() = pt1.Y() - pt3.Y();
    v2 = (pt2 - pt4);
    //	v2.X() = pt2.X() - pt4.X();
    //	v2.Y() = pt2.Y() - pt4.Y();
  } 
  else
  {
    v1 = pt1 - pt2;
    v2 = pt3 - pt4;
    //	v1.X() = pt1.X() - pt2.X();
    //	v1.Y() = pt1.Y() - pt2.Y();
    //	v2.X() = pt3.X() - pt4.X();
    //	v2.Y() = pt3.Y() - pt4.Y();		
  }
/* L. Broglia  
  v3.X(v1.X() + v2.X());
  v3.Y(v1.Y() + v1.Y());
*/
  v3 = v1 + v2 ;
  
  smin =  1.0e8;
  smax = -1.0e8;
  
  G4double norm = sqrt(v3.x() * v3.x() + v3.y() * v3.y());
  if(!norm)
  {
    G4cout << "\nNormal zero!";
    G4cout << "\nLINE & DIR: " << line.x() << " " << line.y() << "  " << dir; 
    G4cout << "\n";
    
    if((abs(line.x())) > kCarTolerance) 
      line.setX(-line.x());
    else
      if((abs(line.y())) > kCarTolerance)
	line.setY(-line.y());
      else
      {
	G4cout << "\n  RETURNING FROm CLIP..";
	smin = 0; smax = 1;
	return;
      }

    G4cout << "\nCHANGED LINE & DIR: " << line.x() << " " 
	   << line.y() << "  " << dir;		
  }
  else
  {
    line.setX( v3.y() / norm);
    line.setY(-v3.x() / norm);
  }

  //    smin =  1.0e8;
  //    smax = -1.0e8;
  //	G4cout << "\n  FINAL LINE & DIR: " << line.X() << " " 
  //           << line.Y() << "  " << dir;	
  
  if( dir == ROW)
  {
    // Create a Convex() hull List 
    for(G4int a = 0; a < col_size; a++)
    {
      ch_ptr = new G4ConvexHull(a/(col_size - 1.0),1.0e8,-1.0e8);
      if(! a) 
      {
	ch_first=ch_ptr;ch_tmp=ch_ptr;
      }
      else ch_tmp->SetNextHull(ch_ptr);
      
      ch_tmp=ch_ptr;
    }
    
    ch_ptr=ch_first;
    register G4double value;
    
    // Loops through the control point mesh and calculates
    // the nvex() hull for the surface.
    
    for( G4int h = 0; h < row_size; h++)
    {
      for(G4int k = 0; k < col_size; k++)
      {
/* L. Broglia
	G4Point2d& coordstmp = (G4Point2d&)ctl_points->get(h,k);  
   	value = - ((coordstmp.X() * line.X() + coordstmp.Y() * line.Y()));
*/
	G4Point3D coordstmp = ctl_points->Get3D(h,k);  
   	value = - ((coordstmp.x() * line.x() + coordstmp.y() * line.y()));

	if( value <= (ch_ptr->GetMin()+kCarTolerance)) ch_ptr->SetMin(value);
	if( value >= (ch_ptr->GetMax()-kCarTolerance)) ch_ptr->SetMax(value);
	    
	ch_ptr=ch_ptr->GetNextHull();
      }
      
      ch_ptr=ch_first;
    }
    
    ch_ptr=ch_first;
    // Finds the points where the nvex() hull intersects
    // with the coordinate .X()is. These points are the
    // minimum and maximum values to where to clip the
    // surface.
    
    for(G4int l = 0; l < col_size - 1; l++)
    {
      ch_tmp=ch_ptr->GetNextHull();
      for(G4int m = l+1; m < col_size; m++)
      {
	register G4double d;
	register G4double param1, param2;
	param1 = ch_ptr->GetParam();
	param2 = ch_tmp->GetParam();
	
	if(ch_tmp->GetMax() - ch_ptr->GetMax())
	{
	  d = Findzero( param1, param2, ch_ptr->GetMax(), ch_tmp->GetMax());
	  if( d <= (smin + kCarTolerance) ) smin = d * .99;
	  if( d >= (smax - kCarTolerance) ) smax = d * .99 + .01;
	}
	
	if(ch_tmp->GetMin() - ch_ptr->GetMin())
	{
	  d = Findzero( param1, param2, ch_ptr->GetMin(), ch_tmp->GetMin());
	  if( d <= (smin + kCarTolerance)) smin = d * .99;
	  if( d >= (smax - kCarTolerance)) smax = d * .99 + .01;
	}
	
	ch_tmp=ch_tmp->GetNextHull();
      }
      
      ch_ptr=ch_ptr->GetNextHull();
    }
    
    ch_ptr=ch_first;

    if (smin <= 0.0)   smin = 0.0;
    if (smax >= 1.0)   smax = 1.0;

    if ( Sign(ch_ptr->GetMin()) != Sign(ch_ptr->GetMax()))  smin = 0.0;
    
    i = Sign(ch_tmp->GetMin()); // ch_tmp points to last nvex()_hull in List
    j = Sign(ch_tmp->GetMax());
    
    if ( abs(i-j) > kCarTolerance ) smax = 1.0;
    //	if ( i != j)  smax = 1.0;
    
  } 
  else // Other G4Vector3D
  {
    for(G4int n = 0; n < row_size; n++)
      {
	ch_ptr = new G4ConvexHull(n/(row_size - 1.0),1.0e8,-1.0e8);
	if(!n) 
	{
	  ch_first=ch_ptr;
	  ch_tmp=ch_ptr;
	}
	else ch_tmp->SetNextHull(ch_ptr);
	
	ch_tmp=ch_ptr;
      }
    
    ch_ptr=ch_first;
    
    for( G4int o = 0; o < col_size; o++)
    {
      for(G4int p = 0; p < row_size; p++)
      {
	register G4double value;

/* L. Broglia
	G4Point2d& coordstmp =(G4Point2d&) ctl_points->get(p,o);  	      
	value = - ((coordstmp.X() * line.X() + coordstmp.Y() * line.Y()));
*/
	G4Point3D coordstmp = ctl_points->Get3D(p,o);  	      
	value = - ((coordstmp.x() * line.x() + coordstmp.y() * line.y()));

	if( value <= (ch_ptr->GetMin()+kCarTolerance)) ch_ptr->SetMin(value);
	if( value >= (ch_ptr->GetMax()-kCarTolerance)) ch_ptr->SetMax(value);
	
	ch_ptr=ch_ptr->GetNextHull();
      }

      ch_ptr=ch_first;
    }
    
    ch_ptr=ch_first;
    ch_tmp=ch_first;
    
    for(G4int q = 0; q < row_size - 1; q++)
    {
      ch_tmp=ch_ptr->GetNextHull();
      for(G4int r = q+1; r < row_size; r++)
      {
	register G4double param1 = ch_ptr->GetParam();
	register G4double param2 = ch_tmp->GetParam();
	register G4double d;
	
	if(ch_tmp->GetMax() - ch_ptr->GetMax())
	{
	  d = Findzero( param1, param2, ch_ptr->GetMax(), ch_tmp->GetMax());
	  if( d <= (smin + kCarTolerance) ) smin = d * .99;
	  if( d >= (smax - kCarTolerance) ) smax = d * .99 + .01;
	}

	if(ch_tmp->GetMin()-ch_ptr->GetMin())
	{
	  d = Findzero( param1, param2, ch_ptr->GetMin(), ch_tmp->GetMin());
	  if( d <= (smin + kCarTolerance) ) smin = d * .99;
	  if( d >= (smax - kCarTolerance) ) smax = d * .99 + .01;
	}
	
	ch_tmp=ch_tmp->GetNextHull();
      }

      ch_ptr=ch_ptr->GetNextHull();
      }
    
    ch_tmp=ch_ptr;
    ch_ptr=ch_first;
	
    if (smin <= 0.0)  smin = 0.0;
    if (smax >= 1.0)  smax = 1.0;
    
    if ( Sign(ch_ptr->GetMin()) != Sign(ch_ptr->GetMax())) smin = 0.0;
    
    i = Sign(ch_tmp->GetMin()); // ch_tmp points to last nvex()_hull in List
    j = Sign(ch_tmp->GetMax());

    //
    if ( (abs(i-j) > kCarTolerance)) smax = 1.0;
  }

  ch_ptr=ch_first;
  while(ch_ptr!=ch_ptr->GetNextHull())
  {
    ch_tmp=ch_ptr;
    ch_ptr=ch_ptr->GetNextHull();
    delete ch_tmp;
  }

  delete ch_ptr;
  
  // Testing...    
  Clips++; 
}


void G4BezierSurface::GetClippedRegionFromSurface()
{
  // Returns the clipped part of the surface. First calculates the
  // length of the new knotvector. Then uses the refinement function to 
  // get the new knotvector and controlmesh.

  //    G4cout << "\nBezier region clipped.";
    
  delete new_knots;
  if ( dir == ROW) 
  {
    new_knots = new G4KnotVector(GetOrder(0) * 2);
    for (register G4int i = 0; i < GetOrder(0); i++) 
    {
      new_knots->PutKnot(i, smin);
      new_knots->PutKnot(i+ GetOrder(0), smax);
    }
  }
  else
  {
    new_knots = new G4KnotVector( GetOrder(1) * 2);
    for ( register G4int i = 0; i <  GetOrder(1); i++) 
    {
      new_knots->PutKnot(i, smin);
      new_knots->PutKnot(i+ GetOrder(1), smax);
    }
  }
} // NURB_REGION_FROM_SURFACE


void G4BezierSurface::RefineSurface()
{
  // Returns the new clipped surface. Calculates the new controlmesh
  // and knotvectorvalues for the surface by using the Oslo-algorithm
  
  delete old_points;
  if (dir == ROW) 
  {
    // Row (u) G4Vector3D 
    ord = GetOrder(0);
    CalcOsloMatrix();
    for(register G4int a=0;a<new_knots->GetSize();a++)
      u_knots->PutKnot(a, new_knots->GetKnot(a));
	
    lower = 0; 
    upper = new_knots->GetSize() - GetOrder(0);
  
    // Copy of the old points.
    old_points = new G4ControlPoints(*ctl_points);
    MapSurface(this);
  }
  else 	
  {	
    ord = GetOrder(1);
    CalcOsloMatrix ();
    for(register G4int a=0;a < new_knots->GetSize();a++)
      v_knots->PutKnot(a, new_knots->GetKnot(a));
	
    // Copy of the old points.
    old_points = new G4ControlPoints(*ctl_points);
    
    // Make new controlpoint matrix,
    register G4int cols = ctl_points->GetCols();
    delete ctl_points;

    ctl_points = new G4ControlPoints(2,(new_knots->GetSize()-
					GetOrder(1)),cols);   
    lower = 0; 
    upper = new_knots->GetSize() - GetOrder(1);
    MapSurface(this);
  }
}// REFINE_SURFACE


void G4BezierSurface::CalcOsloMatrix()
{
  // This algorithm is described in the paper "Making the Oslo-algorithm
  // more efficient" in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86
  // Calculates the oslo-matrix , which is used in mapping the new
  // knotvector- and controlpoint-values.
 
  register G4KnotVector *ah;
  register G4KnotVector *newknots;		     
  register G4int         i;
  register G4int         j;
  register G4int         mu, muprim;
  register G4int         vv, p;
  register G4int         iu, il, ih, n1;		
  register G4int         ahi;	
  register G4double      beta1;
  register G4double      tj;

  ah       = new G4KnotVector(ord*(ord + 1)/2);
  newknots = new G4KnotVector(ord * 2 );
  
  n1 = new_knots->GetSize() - ord;
  mu = 0;		
  
  if(oslo_m!=(G4OsloMatrix*)0)
  {
    G4OsloMatrix* tmp;
    
    //	    while(oslo_m!=oslo_m->next)
    while(oslo_m!=(G4OsloMatrix*)0)	    
    {
      tmp=oslo_m->GetNextNode();delete oslo_m; oslo_m=tmp;
    }
  }
	
  delete oslo_m;
  oslo_m = new G4OsloMatrix();
  
  register G4OsloMatrix* o_ptr = oslo_m;
  
  register G4KnotVector* old_knots;
  if(dir)
    old_knots = v_knots;
  else
    old_knots = u_knots;
  
  for (j = 0; j < n1; j++) 
  {
    if ( j != 0 )
    {
      oslo_m->SetNextNode(new G4OsloMatrix());
      oslo_m = oslo_m->GetNextNode();
    }
    
    while (old_knots->GetKnot(mu + 1) <= new_knots->GetKnot(j))
      mu = mu + 1;		// find the bounding mu 
    
    i = j + 1;
    muprim = mu;
    
    while ((new_knots->GetKnot(i) == old_knots->GetKnot(muprim)) && 
	   i < (j + ord)) 
    {
      i++;
      muprim--;
    }
    
    ih = muprim + 1;
    
    for (vv = 0, p = 1; p < ord; p++) 
    {
      if (new_knots->GetKnot(j + p) == old_knots->GetKnot(ih))
	ih++;
      else
	newknots->PutKnot(++vv - 1,new_knots->GetKnot(j + p));
    }
    
    ahi = AhIndex(0, ord - 1,ord);
    ah->PutKnot(ahi, 1.0);
    
    for (p = 1; p <= vv; p++) 
    {
      beta1 = 0.0;
      tj = newknots->GetKnot(p-1);
      
      if (p - 1 >= muprim) 
      {
	beta1 = AhIndex(p - 1, ord - muprim,ord);
	beta1 = ((tj - old_knots->GetKnot(0)) * beta1) /
	  (old_knots->GetKnot(p + ord - vv) - old_knots->GetKnot(0));
      }

      i  = muprim - p + 1;
      il = Amax (1, i);
      i  = n1 - 1 + vv - p;
      iu = Amin (muprim, i);
      
      for (i = il; i <= iu; i++) 
      {
	register G4double d1, d2;
	register G4double beta;
	
	d1 = tj - old_knots->GetKnot(i);
	d2 = old_knots->GetKnot(i + p + ord - vv - 1) - tj;

	beta = ah->GetKnot(AhIndex(p - 1, i + ord - muprim - 1,ord)) / 
	  (d1 + d2);
				
	
	ah->PutKnot(AhIndex(p, i + ord - muprim - 2,ord), d2 * beta + beta1) ; 
	beta1 = d1 * beta;
      }
      
      ah->PutKnot(AhIndex(p, iu + ord - muprim - 1,ord), beta1);

      if (iu < muprim) 
      {
	register G4double kkk;
	register G4double ahv;
	
	kkk = old_knots->GetKnot(n1 - 1 + ord);
	ahv = AhIndex (p - 1, iu + ord - muprim,ord); 
	ah->PutKnot(AhIndex(p, iu + ord - muprim - 1,ord),
		    beta1 + (kkk - tj) * ahv /
		    (kkk - old_knots->GetKnot(iu + 1)));
      }
    }

    // Remove the oslo matrix List
    G4OsloMatrix* temp_oslo = oslo_m;

/*
      if(oslo_m != (G4OsloMatrix*)0)
      while(oslo_m->next != oslo_m)
      {
      oslo_m = oslo_m->next;
      delete temp_oslo;
      temp_oslo = oslo_m;
      }
      
      // Remove the last
      delete oslo_m;
*/

    while(oslo_m != (G4OsloMatrix*)0)
    {
      oslo_m = oslo_m->GetNextNode();
      delete temp_oslo;
      temp_oslo = oslo_m;
    }
    
    delete oslo_m;
    
    // Create a new oslo matrix    
    oslo_m = new G4OsloMatrix(vv+1, Amax(muprim - vv,0), vv);
    
    for ( i = vv, p = 0; i >= 0; i--)
      oslo_m->GetKnotVector()
            ->PutKnot ( p++, ah->GetKnot(AhIndex (vv, (ord-1) - i,ord)));
    
  }

  delete ah;
  delete newknots;
  oslo_m->SetNextNode(0);
  oslo_m = o_ptr;
}


void G4BezierSurface::MapSurface(G4Surface* tmp)
{
  // This algorithm is described in the paper Making the Oslo-algorithm
  // more efficient in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86
  // Maps the new controlpoints into the new surface.
  
  register G4ControlPoints *c_ptr;
  register G4OsloMatrix *o_ptr;
  register G4ControlPoints* new_pts;
  register G4ControlPoints* old_pts;
  
  new_pts = ctl_points;
  
  // Copy the old points so they can be used in calculating the new ones.
  //  old_pts = new G4ControlPoints(*ctl_points);
  old_pts = old_points;
  register G4int j, 			//	 j loop 
                 i;			//	 oslo loop 

  c_ptr = new_pts;
  register G4int size; // The number of rows or columns, 
                       // depending on processing order

  if(!dir)
    size=new_pts->GetRows();
  else
    size=new_pts->GetCols();

  for(G4int a=0; a<size;a++)
  {
    if ( lower != 0)
      for ( i = 0,  o_ptr = oslo_m; 
	    i < lower; 
	    i++,  o_ptr = o_ptr->GetNextNode());
    else
      o_ptr = oslo_m;
    
    if(!dir)// Direction ROW
    {
      for ( j = lower; j < upper; j++, o_ptr = o_ptr->GetNextNode()) 
      {
	register G4double o_scale;
	register G4int x;
	x=a;

/* L. Broglia	
	register G4Point2d o_pts= (G4Point2d&)old_pts->Get2d(x, o_ptr->GetOffset());
	register G4Point2d tempc= (G4Point2d&)c_ptr->Get2d(j/upper,
							   (j)%upper-lower);
*/
	register G4Point3D o_pts = old_pts->Get3D(x, o_ptr->GetOffset());
	register G4Point3D tempc = c_ptr->Get3D(j/upper, (j)%upper-lower);

	o_scale = o_ptr->GetKnotVector()->GetKnot(0);
	
	tempc.setX(o_pts.x() * o_scale);
	tempc.setY(o_pts.x() * o_scale);
	
	for ( i = 1; i <= o_ptr->GetSize(); i++)
	{
	  o_scale = o_ptr->GetKnotVector()->GetKnot(i);

/* L. Broglia
	  o_pts = (G4Point2d&)old_pts->get(x, i+o_ptr->GetOffset());
	  tempc.X(tempc.X() + o_scale * o_pts.X());
	  tempc.Y(tempc.Y() + o_scale * o_pts.Y());
*/
	  o_pts = old_pts->Get3D(x, i+o_ptr->GetOffset());
	  tempc.setX(tempc.x() + o_scale * o_pts.x());
	  tempc.setY(tempc.y() + o_scale * o_pts.y());

	}
	
	c_ptr->put(a,(j)%upper-lower,tempc);		
      }	
    }
    else  // dir = COL
    {
      for ( j = lower; j < upper; j++, o_ptr = o_ptr->GetNextNode())
      {
	register G4double o_scale;
	register G4int x;
	x=a;

/* L. Broglia	
	register G4Point2d o_pts= (G4Point2d&)old_pts->Get2d(o_ptr->GetOffset(), x);
	register G4Point2d tempc = (G4Point2d&)c_ptr->Get2d((j)%upper-lower,
							    j/upper);
*/
	register G4Point3D o_pts = old_pts->Get3D(o_ptr->GetOffset(), x);
	register G4Point3D tempc = c_ptr->Get3D((j)%upper-lower,j/upper);

	o_scale = o_ptr->GetKnotVector()->GetKnot(0);
	
	tempc.setX(o_pts.x() * o_scale);
	tempc.setY(o_pts.y() * o_scale);
	
	for ( i = 1; i <= o_ptr->GetSize(); i++)
	{
	  o_scale = o_ptr->GetKnotVector()->GetKnot(i);
/* L. Broglia	
	  o_pts= (G4Point2d&)old_pts->get(i+o_ptr->GetOffset(),a);
*/
	  o_pts= old_pts->Get3D(i+o_ptr->GetOffset(),a);			
	  tempc.setX(tempc.x() + o_scale * o_pts.x());
	  tempc.setY(tempc.y() + o_scale * o_pts.y());
	}
	
	c_ptr->put((j)%upper-lower,a,tempc);		
      }
    }
  }
}


void G4BezierSurface::SplitNURBSurface()
{
  // Divides the surface in two parts. Uses the oslo-algorithm to calculate
  // the new knotvectors and controlpoints for  the subsurfaces.
  
  //    G4cout << "\nBezier splitted.";
  
  register G4double value;
  register G4int i;
  register G4int k_index;
  G4BezierSurface *srf1, *srf2;
  G4int nr,nc;
  
  if ( dir == ROW )
  {
    value = u_knots->GetKnot((u_knots->GetSize()-1)/2);
    
    for( i = 0; i < u_knots->GetSize(); i++)
      if( value == u_knots->GetKnot(i) )
      {
	k_index = i;
	break;
      }	

    if ( k_index == 0)
    {
      value = ( value + u_knots->GetKnot(u_knots->GetSize() -1))/2.0;
      k_index = GetOrder(ROW);
    }
	
    new_knots = u_knots->MultiplyKnotVector(GetOrder(ROW), value);
    
    ord = GetOrder(ROW);
    CalcOsloMatrix();
	
    srf1 = new G4BezierSurface(*this);
    //	srf1->dir=ROW;
    srf1->dir=COL;	
    
    new_knots->ExtractKnotVector(srf1->u_knots, k_index +
				 srf1->GetOrder(ROW),0); 

    nr= srf1->v_knots->GetSize() - srf1->GetOrder(COL);
    nc= srf1->u_knots->GetSize() - srf1->GetOrder(ROW);
    delete srf1->ctl_points;
    
    srf1->ctl_points= new G4ControlPoints(2, nr, nc);
    srf2 = new G4BezierSurface(*this);

    //	srf2->dir = ROW;
    srf2->dir = COL;	

    new_knots->ExtractKnotVector(srf2->u_knots, 
				 new_knots->GetSize(), k_index); 
    
    nr= srf2->v_knots->GetSize() - srf2->GetOrder(COL);
    nc= srf2->u_knots->GetSize() - srf2->GetOrder(ROW);
    
    delete  srf2->ctl_points;
    srf2->ctl_points = new G4ControlPoints(2, nr, nc);
    
    lower = 0;
    upper = k_index;
    MapSurface(srf1);
    
    lower = k_index;
    upper = new_knots->GetSize() - srf2->GetOrder(ROW);
    MapSurface(srf2);
  }
  else // G4Vector3D = col
  {
    value = v_knots->GetKnot((v_knots->GetSize() -1)/2);
    
    for( i = 0; i < v_knots->GetSize(); i++)
      if( value == v_knots->GetKnot(i))
      {
	k_index = i;
	break;
      }
    if ( k_index == 0)
    {
      value = ( value + v_knots->GetKnot(v_knots->GetSize() -1))/2.0;
      k_index = GetOrder(COL);
    }
    
    new_knots = v_knots->MultiplyKnotVector( GetOrder(COL), value );
    ord = GetOrder(COL);
    
    CalcOsloMatrix();
    
    srf1 = new G4BezierSurface(*this);
    //	srf1->dir = COL;
    srf1->dir = ROW;
    
    new_knots->ExtractKnotVector(srf1->v_knots, 
				 k_index + srf1->GetOrder(COL), 0);
	
    nr = srf1->v_knots->GetSize() - srf1->GetOrder(COL);
    nc = srf1->u_knots->GetSize() - srf1->GetOrder(ROW);
    
    delete srf1->ctl_points;
    srf1->ctl_points = new G4ControlPoints(2, nr, nc);
    
    srf2 = new G4BezierSurface(*this);
    //	srf2->dir = COL;
    srf2->dir = ROW;
    
    new_knots->ExtractKnotVector(srf2->v_knots, new_knots->GetSize(), k_index);
	
    nr = srf2->v_knots->GetSize() - srf2->GetOrder(COL);
    nc = srf2->u_knots->GetSize() - srf2->GetOrder(ROW);
    
    delete srf2->ctl_points;
    srf2->ctl_points = new G4ControlPoints(2,nr, nc);
    
    lower = 0;
    upper = k_index; 
    MapSurface(srf1);

    //	next->oslo_m = oslo_m;
    lower = k_index;
    upper = new_knots->GetSize() - srf2->GetOrder(COL);
    MapSurface(srf2);
  }
  
  bezier_list->AddSurface(srf1);
  bezier_list->AddSurface(srf2);
  delete new_knots;
  
  // Testing
  Splits++;  
}
