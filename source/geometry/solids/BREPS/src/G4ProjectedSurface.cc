// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProjectedSurface.cc,v 1.1 1999-01-07 16:07:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ProjectedSurface.hh"
 
int G4ProjectedSurface::Splits=0;

G4ProjectedSurface::~G4ProjectedSurface()
{
  delete u_knots;
  delete v_knots;
  delete ctl_points;
  
  G4OsloMatrix* temp_oslo;
  if(oslo_m!=(G4OsloMatrix*)0)
  {
    while(oslo_m->next != oslo_m)
    {
      temp_oslo = oslo_m;
      oslo_m    = oslo_m->next;
      
      delete temp_oslo;
    }
    
    delete oslo_m;
  }
  
  delete bbox;
}


G4ProjectedSurface::G4ProjectedSurface()
{
  distance = 0;
  oslo_m   =(G4OsloMatrix*)0;
}


G4ProjectedSurface::G4ProjectedSurface(const G4ProjectedSurface &tmp) 
{
  distance = tmp.distance;
  oslo_m   =(G4OsloMatrix*)0;  
   
  //    next=this;
  
  order[0] = tmp.order[0];
  order[1] = tmp.order[1];
  dir = tmp.dir;
  
  u_knots = new G4KnotVector(*tmp.u_knots);
  v_knots = new G4KnotVector(*tmp.v_knots);
  
  ctl_points = new G4ControlPoints(*tmp.ctl_points); 
  bbox = new G4BoundingBox3D();
}


void  G4ProjectedSurface::CopySurface()
  // Copies the projected surface into a bezier surface
  // and adds it to the List of bezier surfaces.
{
  G4BezierSurface *bez = new G4BezierSurface();
  
  bez->Distance(distance);
  bez->PutOrder(0, order[0]);
  bez->PutOrder(1, order[1]);
  bez->Dir(dir);
  bez->u_knots = new G4KnotVector(*u_knots);
  bez->v_knots = new G4KnotVector(*v_knots);
  bez->ctl_points = new G4ControlPoints(*ctl_points);
  
  bezier_list->AddSurface(bez);
}


void G4ProjectedSurface::CalcBBox()
{
  // Finds the bounds of the 2D-projected nurb iow
  // calculates the bounds for a bounding rectangle
  // to the surface. The bounding rectangle is used
  // for a preliminary check of intersection.
  
  // Loop to search the whole control point mesh
  // for the minimum and maximum values for x and y.
  G4double box_minx,box_miny,box_maxx,box_maxy;
  box_minx = INFINITY;
  box_miny = INFINITY;
  box_maxx  = -INFINITY;
  box_maxy  = -INFINITY;
    
  G4double bminx,bminy,bmaxx,bmaxy,tmpx,tmpy;
  bminx = box_minx; bminy = box_miny;
  bmaxx = box_maxx; bmaxy = box_maxy;	    

  for(register int a = ctl_points->GetRows()-1; a>=0;a--)
    for(register int b = ctl_points->GetCols()-1; b>=0;b--)
    {	    
/* L. Broglia
      G4Point2d& tmp = (G4Point2d&)ctl_points->get(a,b);
*/
      G4Point3D tmp = ctl_points->Get3D(a,b);

      tmpx = tmp.x(); tmpy = tmp.y();
      if(bminx > tmpx) box_minx=tmpx;
      if(bmaxx < tmpx) box_maxx=tmpx;	
      if(bminy > tmpy) box_miny=tmpy;	
      if(bmaxy < tmpy) box_maxy=tmpy;
    }
    
  G4Point3D box_min(box_minx,box_miny);
  G4Point3D box_max(box_maxx,box_maxy);

  delete bbox;
  bbox = new G4BoundingBox3D(box_min, box_max);
}


void G4ProjectedSurface::ConvertToBezier(G4SurfaceList& proj_list,
					 G4SurfaceList& bez_list)
{
  projected_list = &proj_list;
  bezier_list    = &bez_list;
  
  // Check wether the surface is a bezier surface by checking
  // if internal knots exist.
  if(CheckBezier())
  {
    // Make it a G4BezierSurface -object and add it to the bezier
    // surface List	
    CopySurface();  
	
    // Retrieve a pointer to the newly added surface iow the
    // last in the List
    G4BezierSurface* bez_ptr = (G4BezierSurface*)bezier_list->GetLastSurface();
    
    // Do the first clip to the bezier.
    bez_ptr->ClipSurface();
    G4double dMin =  bez_ptr->SMin();
    G4double dMax =  bez_ptr->SMax();
    G4double dMaxMinusdMin = dMax - dMin;
    
    if(( dMaxMinusdMin > kCarTolerance ))
    {
      if( dMaxMinusdMin > 0.8 )
      {
	// The clipping routine selected a larger Area than one
	// knot interval which indicates that we have a case of
	// multiple intersections. The projected surface has to
	// be split again in order to separate the intersections
	// to different surfaces.

	// Check tolerance of clipping
	//	    G4cout << "\nClip Area too big -> Split";
	dir = bez_ptr->dir;
	bezier_list->RemoveSurface(bez_ptr);
	
	SplitNURBSurface();
	return;
	//}
      }
      else
	if( dMin > 0.0 || dMax < 0.0 )
	{
	  // The ray intersects with the bounding box
	  // but not with the surface itself.  
	  //		G4cout << "\nConvex hull missed.";
	  bezier_list->RemoveSurface(bez_ptr);
	  return;
	}
    }
    else
      if(dMaxMinusdMin < kCarTolerance && dMaxMinusdMin > -kCarTolerance)
      {
	bezier_list->RemoveSurface(bez_ptr);
	return;
      }

    bez_ptr->LocalizeClipValues();
    bez_ptr->SetValues();	

    // Other G4ThreeVec clipping and testing.
    bez_ptr->ChangeDir();//bez->dir = !bez_ptr->dir;
    bez_ptr->ClipSurface();
    //	G4cout<<"\nSMIN: " <<  bez_ptr->smin << "  SMAX: " 
    //        <<  bez_ptr->smax << " DIR: " << bez_ptr->dir;
	   
    dMin = bez_ptr->SMin();
    dMax = bez_ptr->SMax();
    dMaxMinusdMin = dMax-dMin;	
    
    if((dMaxMinusdMin > kCarTolerance ))// ||
      //	   (dMaxMinusdMin < -kCarTolerance))
    {
      if( (dMaxMinusdMin) > 0.8 )
      {
	//	    G4cout << "\nClip Area too big -> Split";	    
	dir = bez_ptr->dir;//1.2 klo 18.30
	
	//	    dir=!dir;
	bezier_list->RemoveSurface(bez_ptr);
	SplitNURBSurface();
	return;
		 //}
      }
      else
	if( dMin > 1.0 || dMax < 0.0 )
	{
	  //		G4cout << "\nConvex hull missed.";
	  bezier_list->RemoveSurface(bez_ptr);
	  return;
	}
    }
    else
      if(dMaxMinusdMin < kCarTolerance && dMaxMinusdMin > -kCarTolerance)
      {
	bezier_list->RemoveSurface(bez_ptr);
	return;
      }
    
    bez_ptr->LocalizeClipValues();	
    bez_ptr->SetValues();
    bez_ptr->CalcAverage();
  }
  else
  {
    // Split the surface into two new surfaces. The G4ThreeVec
    // is set in the CheckBezier function. 
    //	G4cout << "\nNot a bezier surface -> Split";	    	
    SplitNURBSurface();
  }
}


int G4ProjectedSurface::CheckBezier()
{
  // Checks if the surface is a bezier surface by
  // checking wether internal knots exist. If no internal
  // knots exist the quantity of knots is 2*order of the
  // surface. Returns 1 if the surface 
  // is a bezier.
  
  if( u_knots->GetSize() > (2.0 * GetOrder(ROW)))
	{dir=0;return 0;}
  
  if( v_knots->GetSize() > (2.0 * GetOrder(COL)))
    {dir=1;return 0;}
  
  return 1;
}


void G4ProjectedSurface::SplitNURBSurface()
{
  // Divides the surface in two parts. Uses the oslo-algorithm to calculate
  // the new knotvectors and controlpoints for  the subsurfaces.

  //    G4cout << "\nProjected splitted.";
    
  register G4double value;
  register int i;
  register int k_index;
  register G4ProjectedSurface *srf1, *srf2;
  register int nr,nc;
    
  if ( dir == ROW )
  {
    value = u_knots->GetKnot((u_knots->GetSize()-1)/2);
    
    for( i = 0; i < u_knots->GetSize(); i++)
      if( (abs(value - u_knots->GetKnot(i))) < kCarTolerance )
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
    
    srf1 = new G4ProjectedSurface(*this);
	
    //srf1->dir=ROW;
    srf1->dir=COL;	

    new_knots->ExtractKnotVector(srf1->u_knots, 
				 k_index + srf1->GetOrder(ROW),0); 

    nr= srf1->v_knots->GetSize() - srf1->GetOrder(COL);
    nc= srf1->u_knots->GetSize() - srf1->GetOrder(ROW);
    
    delete srf1->ctl_points;
    srf1->ctl_points= new G4ControlPoints(2, nr, nc);
    
    srf2 = new G4ProjectedSurface(*this);
    
    //srf2->dir = ROW;
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
  else // G4ThreeVec = col
  {
    value = v_knots->GetKnot((v_knots->GetSize() -1)/2);
    
    for( i = 0; i < v_knots->GetSize(); i++)
      if( (abs(value - v_knots->GetKnot(i))) < kCarTolerance )
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
	
    srf1 = new G4ProjectedSurface(*this);

    //srf1->dir = COL;
    srf1->dir = ROW;
    
    new_knots->ExtractKnotVector(srf1->v_knots, 
				 k_index + srf1->GetOrder(COL), 0);
    
    nr = srf1->v_knots->GetSize() - srf1->GetOrder(COL);
    nc = srf1->u_knots->GetSize() - srf1->GetOrder(ROW);
	
    delete srf1->ctl_points;
    srf1->ctl_points = new G4ControlPoints(2, nr, nc);
    
    srf2 = new G4ProjectedSurface(*this);

    //srf2->dir = COL;
    srf2->dir = ROW;
    
    new_knots->ExtractKnotVector(srf2->v_knots, new_knots->GetSize(), k_index);
	
    nr = srf2->v_knots->GetSize() - srf2->GetOrder(COL);
    nc = srf2->u_knots->GetSize() - srf2->GetOrder(ROW);
    
    delete srf2->ctl_points;
    srf2->ctl_points = new G4ControlPoints(2,nr, nc);
    
    lower = 0;
    upper = k_index; 
    MapSurface(srf1);
    
    lower = k_index;
    upper = new_knots->GetSize() - srf2->GetOrder(COL);
    MapSurface(srf2);  
  }
  
  // Check that surfaces are ok.
  int col_size = srf1->ctl_points->GetCols();
  int row_size = srf1->ctl_points->GetRows();

/* L. Broglia  
  // get three cornerpoints of the controlpoint mesh.
  G4Point2d pt1 = srf1->ctl_points->get(0,0);    
  G4Point2d pt2 =  srf1->ctl_points->get(0,col_size-1);    
  G4Point2d pt3 =  srf1->ctl_points->get(row_size-1,0);    
  
  // Calc distance between points
  G4double pointDist1 = pt1.Distance(pt2);
  G4double pointDist2 = pt1.Distance(pt3);
*/

  // get three cornerpoints of the controlpoint mesh.
  G4Point3D pt1 =  srf1->ctl_points->Get3D(0,0);    
  G4Point3D pt2 =  srf1->ctl_points->Get3D(0,col_size-1);    
  G4Point3D pt3 =  srf1->ctl_points->Get3D(row_size-1,0);    
  
  // Calc distance squared between points
  G4double pointDist1 = pt1.distance2(pt2);
  G4double pointDist2 = pt1.distance2(pt3);


  // Add surfaces to List of projected surfaces    
  if(pointDist1 > kCarTolerance && pointDist2 > kCarTolerance)
    projected_list->AddSurface(srf1);
  else
    delete srf1;

  col_size = srf2->ctl_points->GetCols();
  row_size = srf2->ctl_points->GetRows();

/* L. Broglia      
  // get three cornerpoints of the controlpoint mesh.
  pt1 = srf2->ctl_points->get(0,0);    
  pt2 =  srf2->ctl_points->get(0,col_size-1);    
  pt3 =  srf2->ctl_points->get(row_size-1,0);    
  
  // Calc distance between points
  pointDist1 = pt1.Distance(pt2);
  pointDist2 = pt1.Distance(pt3);
*/

  // get three cornerpoints of the controlpoint mesh.
  pt1 =  srf2->ctl_points->Get3D(0,0);    
  pt2 =  srf2->ctl_points->Get3D(0,col_size-1);    
  pt3 =  srf2->ctl_points->Get3D(row_size-1,0);    
  
  // Calc distance squared between points
  pointDist1 = pt1.distance2(pt2);
  pointDist2 = pt1.distance2(pt3);

  // Add surfaces to List of projected surfaces    
  if(pointDist1 > kCarTolerance && pointDist2 > kCarTolerance)
    projected_list->AddSurface(srf2);
  else
    delete srf2;
  
  delete new_knots;
  
  Splits++;   
}

void G4ProjectedSurface::CalcOsloMatrix()
{
  // This algorithm is described in the paper "Making the Oslo-algorithm
  // more efficient" in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86
  // Calculates the oslo-matrix , which is used in mapping the new
  // knotvector- and controlpoint-values.
 
  register G4KnotVector *ah;
  static G4KnotVector *newknots;		     
  register int     i;
  register int      j;
  register int      mu, muprim;
  register int      v, p;
  register int      iu, il, ih, n1;		
  register int      ahi;	
  register G4double beta1;
  register G4double tj;
	
  ah = new G4KnotVector(ord*(ord + 1)/2);
	
  newknots = new G4KnotVector(ord * 2 );

  n1 = new_knots->GetSize() - ord;
  mu = 0;		
  
  if(oslo_m!=(G4OsloMatrix*)0)
  {
    G4OsloMatrix* tmp;
    while(oslo_m!=oslo_m->next)
    {
      tmp=oslo_m->next;
      delete oslo_m; 
      oslo_m=tmp;
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
      oslo_m->next = new G4OsloMatrix();
      oslo_m = oslo_m->next;
    }

    //while (old_knots->GetKnot(mu + 1) <= new_knots->GetKnot(j))
    while ( (new_knots->GetKnot(j) - old_knots->GetKnot(mu + 1)) > 
	    kCarTolerance                                            )
      mu = mu + 1;		// find the bounding mu 
    
    i = j + 1;
    muprim = mu;
    
    while ( ((abs(new_knots->GetKnot(i) - old_knots->GetKnot(muprim))) < 
	     kCarTolerance) && i < (j + ord)                             ) 
    {
      i++;
      muprim--;
    }

    ih = muprim + 1;
    
    for (v = 0, p = 1; p < ord; p++) 
    {
      // if (new_knots->GetKnot(j + p) == old_knots->GetKnot(ih))
      if ( (abs((new_knots->GetKnot(j + p)) - (old_knots->GetKnot(ih)))) < 
	   kCarTolerance                                                    )
	ih++;
      else
	newknots->PutKnot(++v - 1,new_knots->GetKnot(j + p));
    }

    ahi = AhIndex(0, ord - 1,ord);
    ah->PutKnot(ahi, 1.0);
    
    for (p = 1; p <= v; p++) 
    {
      beta1 = 0.0;
      tj = newknots->GetKnot(p-1);
      
      if (p - 1 >= muprim) 
      {
	beta1 = AhIndex(p - 1, ord - muprim,ord);
	beta1 = ((tj - old_knots->GetKnot(0)) * beta1) /
	  (old_knots->GetKnot(p + ord - v) - old_knots->GetKnot(0));
      }
	
      i  = muprim - p + 1;
      il = Amax (1, i);
      i  = n1 - 1 + v - p;
      iu = Amin (muprim, i);
      
      for (i = il; i <= iu; i++) 
      {
	register G4double d1, d2;
	register G4double beta;
	
	d1 = tj - old_knots->GetKnot(i);
	d2 = old_knots->GetKnot(i + p + ord - v - 1) - tj;
	
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
    
    delete oslo_m->o_vec;
    oslo_m->o_vec = new G4KnotVector(v+1);
    oslo_m->offset = Amax(muprim - v, 0);
    oslo_m->osize = v;
    
    for ( i = v, p = 0; i >= 0; i--)
      oslo_m->o_vec->PutKnot ( p++, ah->GetKnot(AhIndex (v, (ord-1) - i,ord)));
    
  }
  
  delete ah;
  delete newknots;
  oslo_m->next = oslo_m;
  oslo_m = o_ptr;
}














void  G4ProjectedSurface::MapSurface(G4ProjectedSurface* srf)
{
  // This algorithm is described in the paper "Making the Oslo-algorithm
  // more efficient" in SIAM J.NUMER.ANAL. Vol.23, No. 3, June '86
  // Maps the new controlpoints into the new surface.

  register G4ControlPoints *c_ptr;
  register G4OsloMatrix    *o_ptr;
  register G4ControlPoints* new_pts;
  register G4ControlPoints* old_pts;
  
  new_pts = srf->ctl_points;
  
  // Copy the old points so they can be used in calculating the new ones.
  // In this version, where the splitted surfaces are given
  // as parameters the copying is not necessary.
  
  old_pts = new G4ControlPoints(*ctl_points);
  register int	j,    //  j loop 
                i;    //  oslo loop 
  c_ptr = new_pts;
  
  register int size; // The number of rows or columns, 
                     // depending on processing order

  if(!dir)
    size=new_pts->GetRows();
  else
    size=new_pts->GetCols();
  
  for( register int a=0; a<size;a++)
  {
    if ( lower != 0)
      for ( i = 0,  o_ptr = oslo_m; i < lower; i++,  o_ptr = o_ptr->next);
    else
      o_ptr = oslo_m;
    
    if(!dir)// Direction ROW
    {
      for ( j = lower; j < upper; j++, o_ptr = o_ptr->next) 
      {
	register G4double o_scale;
	register int x;
	x=a;

/* L. Broglia	
	register G4Point2d o_pts = (G4Point2d&)old_pts->get(x,o_ptr->offset);
	register G4Point2d tempc = (G4Point2d&)c_ptr->get(j/upper,
							  (j)%upper-lower);
*/
	register G4Point3D o_pts = old_pts->Get3D(x, o_ptr->offset);
	register G4Point3D tempc = c_ptr->Get3D(j/upper, (j)%upper-lower);
	o_scale = o_ptr->o_vec->GetKnot(0);

	tempc.setX(o_pts.x() * o_scale);
	tempc.setY(o_pts.y() * o_scale);

	for ( i = 1; i <= o_ptr->osize; i++) 
	{
	  o_scale = o_ptr->o_vec->GetKnot(i);

/* L. Broglia	  
	  o_pts = (G4Point2d&)old_pts->get(x,i+o_ptr->offset);
	  tempc.X(tempc.X() + o_scale * o_pts.X());
	  tempc.Y(tempc.Y() + o_scale * o_pts.Y());
*/

	  o_pts = old_pts->Get3D(x,i+o_ptr->offset);
	  tempc.setX(tempc.x() + o_scale * o_pts.x());
	  tempc.setY(tempc.y() + o_scale * o_pts.y());
	}
	
	c_ptr->put(a,(j)%upper-lower,tempc);		
      }
    }
    else  // dir = COL
    {
      for ( j = lower; j < upper; j++, o_ptr = o_ptr->next) 
      {
	register G4double o_scale;
	register int x;
	x=a;

/* L.Broglia
	register G4Point2d o_pts = (G4Point2d&)old_pts->get(o_ptr->offset,x);
	register G4Point2d tempc = (G4Point2d&)c_ptr->get((j)%upper-lower,
							  j/upper);
*/
	register G4Point3D o_pts = old_pts->Get3D(o_ptr->offset,x);
	register G4Point3D tempc = c_ptr->Get3D((j)%upper-lower, j/upper);
		
	o_scale = o_ptr->o_vec->GetKnot(0);

	tempc.setX(o_pts.x() * o_scale);
	tempc.setY(o_pts.y() * o_scale);

	for ( i = 1; i <= o_ptr->osize; i++) 
	{
	  o_scale = o_ptr->o_vec->GetKnot(i);
	  o_pts= old_pts->Get3D(i+o_ptr->offset,a);
	  
	  tempc.setX(tempc.x() + o_scale * o_pts.x());
	  tempc.setY(tempc.y() + o_scale * o_pts.y());
	}
	
	c_ptr->put((j)%upper-lower,a,tempc);
      }
    }
  }
  
  delete old_pts;
}













