// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ToroidalSurface.cc,v 1.4 2000-11-08 14:22:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ToroidalSurface.cc
//
// ----------------------------------------------------------------------

#include "G4ToroidalSurface.hh"


G4ToroidalSurface::G4ToroidalSurface()
 : EQN_EPS(1e-9)
{
}

G4ToroidalSurface::G4ToroidalSurface(const G4Vector3D& Location,
				     const G4Vector3D& Ax,
				     const G4Vector3D& Dir,
				     G4double MinRad,
				     G4double MaxRad)
  : EQN_EPS(1e-9)
{   
  Placement.Init(Dir, Ax, Location);

  MinRadius = MinRad;
  MaxRadius = MaxRad;
  TransMatrix= new G4PointMatrix(4,4);
}


G4ToroidalSurface::~G4ToroidalSurface()
{
}


void G4ToroidalSurface::CalcBBox()
{
  // L. Broglia
  // G4Point3D Origin = Placement.GetSrfPoint();
  G4Point3D Origin = Placement.GetLocation();

  G4Point3D Min(Origin.x()-MaxRadius,
		Origin.y()-MaxRadius,
		Origin.z()-MaxRadius);
  G4Point3D Max(Origin.x()+MaxRadius,
		Origin.y()+MaxRadius,
		Origin.z()+MaxRadius);
 
  bbox = new G4BoundingBox3D(Min,Max);
}


G4Vector3D G4ToroidalSurface::SurfaceNormal(const G4Point3D& Pt) const 
{
  return G4Vector3D(0,0,0);
}


G4double G4ToroidalSurface::ClosestDistanceToPoint(const G4Point3D &Pt)
{
  // L. Broglia
  // G4Point3D Origin = Placement.GetSrfPoint();
  G4Point3D Origin = Placement.GetLocation();

  G4double  Dist   = Pt.distance(Origin);

  return ((Dist - MaxRadius)*(Dist - MaxRadius));
}


G4int G4ToroidalSurface::Intersect(const G4Ray& Ray)
{
  // ----	inttor - Intersect a ray with a torus. ------------------------
  //	from GraphicsGems II by 
  
  //	Description:							 
  //	    Inttor determines the intersection of a ray with a torus.	 
  //									 
  //	On entry:							 
  //	    raybase = The coordinate defining the base of the		 
  //		      intersecting ray.					 
  //	    raycos  = The G4Vector3D cosines of the above ray.		 
  //	    center  = The center location of the torus.			 
  //	    radius  = The major radius of the torus.			 
  //	    rplane  = The minor radius in the G4Plane of the torus.	 
  //	    rnorm   = The minor radius Normal to the G4Plane of the torus. 
  //	    tran    = A 4x4 transformation matrix that will position	 
  //		      the torus at the origin and orient it such that	 
  //		      the G4Plane of the torus lyes in the x-z G4Plane.	 
  //									 
  //	On return:							 
  //	    nhits   = The number of intersections the ray makes with	 
  //		      the torus.					 
  //	    rhits   = The entering/leaving distances of the		 
  //		      intersections.					 
  //									 
  //	Returns:  True if the ray intersects the torus.			 
  //									 
  // --------------------------------------------------------------------
	   
  // Variables. Should be optimized later...
  G4Point3D  Base = Ray.GetStart();  // Base of the intersection ray
  G4Vector3D DCos = Ray.GetDir();    // Direction cosines of the ray
  G4int	     nhits=0;		     // Number of intersections
  G4double   rhits[4];		     // Intersection distances
  G4double   hits[4];		     // Ordered intersection distances
  G4double   rho, a0, b0;	     // Related constants		
  G4double   f, l, t, g, q, m, u;    // Ray dependent terms		
  G4double   C[5];		     // Quartic coefficients	      
	
  //	Transform the intersection ray					
  
  
  //	MultiplyPointByMatrix  (Base);  // Matriisi puuttuu viela!
  //	MultiplyVectorByMatrix (DCos);
  
  //	Compute constants related to the torus.	
  G4double rnorm = MaxRadius - MinRadius; // ei tietoa onko oikein...
  rho = MinRadius*MinRadius / (rnorm*rnorm);
  a0  = 4. * MaxRadius*MaxRadius;
  b0  = MaxRadius*MaxRadius - MinRadius*MinRadius;
  
  //	Compute ray dependent terms.				       
  f = 1. - DCos.y()*DCos.y();
  l = 2. * (Base.x()*DCos.x() + Base.z()*DCos.z());
  t = Base.x()*Base.x() + Base.z()*Base.z();
  g = f + rho * DCos.y()*DCos.y();
  q = a0 / (g*g);
  m = (l + 2.*rho*DCos.y()*Base.y()) / g;
  u = (t +    rho*Base.y()*Base.y() + b0) / g;
	
  //	Compute the coefficients of the quartic.			
  
  C[4] = 1.0;
  C[3] = 2. * m;
  C[2] = m*m + 2.*u - q*f;
  C[1] = 2.*m*u - q*l;
  C[0] = u*u - q*t;
	
  //	Use quartic root solver found in "Graphics Gems" by Jochen Schwarze.
  nhits = SolveQuartic (C,rhits);
  
  //	SolveQuartic returns root pairs in reversed order.		
  m = rhits[0]; u = rhits[1]; rhits[0] = u; rhits[1] = m;
  m = rhits[2]; u = rhits[3]; rhits[2] = u; rhits[3] = m;

  //  	return (*nhits != 0);
  
  if(nhits != 0)
  {
    // Convert Hit distances to intersection points
    /*
      G4Point3D** IntersectionPoints = new G4Point3D*[nhits];
      for(G4int a=0;a<nhits;a++)
      {
      G4double Dist = rhits[a];
      IntersectionPoints[a] = new G4Point3D((Base - Dist * DCos)); 
      }
      // Check wether any of the hits are on the actual surface
      // Start with checking for the intersections that are Inside the bbox
      
      G4Point3D* Hit;
      G4int InsideBox[2]; // Max 2 intersections on the surface
      G4int Counter=0;
    */

    G4Point3D BoxMin = bbox->GetBoxMin();
    G4Point3D BoxMax = bbox->GetBoxMax();
    G4Point3D Hit;
    G4int       c1     = 0;
    G4int       c2;
    G4double  tempVec[4];
    
    for(G4int a=0;a<nhits;a++) 
    {
      while ( (c1 < 4) && (hits[c1] <= rhits[a]) )
      {
	tempVec[c1]=hits[c1]; 
	c1++;
      }
	
      for(c2=c1+1;c2<4;c2++) 
	tempVec[c2]=hits[c2-1];
	
      if(c1<4) 
      {
	tempVec[c1]=rhits[a];
	
	for(c2=0;c2<4;c2++)
	  hits[c2]=tempVec[c2];
      }
    }
    
    for(G4int b=0;b<nhits;b++)
    {
      //		Hit = IntersectionPoints[b]; 
      if(hits[b] >=kCarTolerance*0.5)
      {
	Hit = Base + (hits[b]*DCos);
	//	      InsideBox[Counter]=b;
	if( (Hit.x() > BoxMin.x()) &&
	    (Hit.x() < BoxMax.x()) &&
	    (Hit.y() > BoxMin.y()) &&
	    (Hit.y() < BoxMax.y()) &&		
	    (Hit.z() > BoxMin.z()) &&
	    (Hit.z() < BoxMax.z())    )
	{
	  closest_hit = Hit;
	  distance =  hits[b]*hits[b];
	  return 1;
	}
	
	//		    Counter++;
      }
    }
    
    // If two Inside bbox, find closest 
    G4int Closest=0;
    
    //	    if(Counter>1)
    //		if(rhits[InsideBox[0]] > rhits[InsideBox[1]])
    //		    Closest=1;
    
    // Project polygon and do point in polygon
    // Projection also for curves etc.
    // Should probably be implemented in the curve class. 
    G4Plane Plane1 = Ray.GetPlane(1);
    G4Plane Plane2 = Ray.GetPlane(2);
    
    // Point in polygon
    return 1;
  }
  return 0;
}


G4int G4ToroidalSurface::SolveQuartic(G4double c[], G4double s[]  )
{
  // From Graphics Gems I by Jochen Schwartz
  
  G4double  coeffs[ 4 ];
  G4double  z, u, v, sub;
  G4double  A, B, C, D;
  G4double  sq_A, p, q, r;
  G4int     i, num;
  
    // Normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 

  A = c[ 3 ] / c[ 4 ];
  B = c[ 2 ] / c[ 4 ];
  C = c[ 1 ] / c[ 4 ];
  D = c[ 0 ] / c[ 4 ];

  //  substitute x = y - A/4 to eliminate cubic term:
  // x^4 + px^2 + qx + r = 0 
  
  sq_A = A * A;
  p    = - 3.0/8 * sq_A + B;
  q    = 1.0/8 * sq_A * A - 1.0/2 * A * B + C;
  r    = - 3.0/256*sq_A*sq_A + 1.0/16*sq_A*B - 1.0/4*A*C + D;
  
  if (IsZero(r))
  {
    // no absolute term: y(y^3 + py + q) = 0 
    
    coeffs[ 0 ] = q;
    coeffs[ 1 ] = p;
    coeffs[ 2 ] = 0;
    coeffs[ 3 ] = 1;
    
    num = SolveCubic(coeffs, s);
    
    s[ num++ ] = 0;
  }
  else
  {
    // solve the resolvent cubic ... 
    coeffs[ 0 ] = 1.0/2 * r * p - 1.0/8 * q * q;
    coeffs[ 1 ] = - r;
    coeffs[ 2 ] = - 1.0/2 * p;
    coeffs[ 3 ] = 1;
    
    (void) SolveCubic(coeffs, s);
    
    // ... and take the one real solution ... 
    z = s[ 0 ];

    // ... to Build two quadric equations 
    u = z * z - r;
    v = 2 * z - p;
    
    if (IsZero(u))
      u = 0;
    else if (u > 0)
      u = sqrt(u);
    else
      return 0;

    if (IsZero(v))
      v = 0;
    else if (v > 0)
      v = sqrt(v);
    else
      return 0;

    coeffs[ 0 ] = z - u;
    coeffs[ 1 ] = q < 0 ? -v : v;
    coeffs[ 2 ] = 1;
    
    num = SolveQuadric(coeffs, s);
    
    coeffs[ 0 ]= z + u;
    coeffs[ 1 ] = q < 0 ? v : -v;
    coeffs[ 2 ] = 1;
    
    num += SolveQuadric(coeffs, s + num);
  }
  
  // resubstitute 
  
  sub = 1.0/4 * A;
  
  for (i = 0; i < num; ++i)
    s[ i ] -= sub;
  
  return num;
}


G4int G4ToroidalSurface::SolveCubic(G4double c[], G4double s[]  )
{
  // From Graphics Gems I bu Jochen Schwartz
  G4int     i, num;
  G4double  sub;
  G4double  A, B, C;
  G4double  sq_A, p, q;
  G4double  cb_p, D;

  // Normal form: x^3 + Ax^2 + Bx + C = 0 
  A = c[ 2 ] / c[ 3 ];
  B = c[ 1 ] / c[ 3 ];
  C = c[ 0 ] / c[ 3 ];
  
  //  substitute x = y - A/3 to eliminate quadric term:
  //	x^3 +px + q = 0 
  sq_A = A * A;
  p = 1.0/3 * (- 1.0/3 * sq_A + B);
  q = 1.0/2 * (2.0/27 * A * sq_A - 1.0/3 * A * B + C);
  
  // use Cardano's formula 
  cb_p = p * p * p;
  D = q * q + cb_p;
  
  if (IsZero(D))
  {
    if (IsZero(q)) // one triple solution 
    {
      s[ 0 ] = 0;
      num = 1;
    }
    else // one single and one G4double solution 
    {
      G4double u = cbrt(-q);
      s[ 0 ] = 2 * u;
      s[ 1 ] = - u;
      num = 2;
    }
  }
  else if (D < 0) // Casus irreducibilis: three real solutions
  {
    G4double phi = 1.0/3 * acos(-q / sqrt(-cb_p));
    G4double t = 2 * sqrt(-p);
    
    s[ 0 ] =   t * cos(phi);
    s[ 1 ] = - t * cos(phi + M_PI / 3);
    s[ 2 ] = - t * cos(phi - M_PI / 3);
    num = 3;
  }
  else // one real solution 
  {
    G4double sqrt_D = sqrt(D);
    G4double u = cbrt(sqrt_D - q);
    G4double v = - cbrt(sqrt_D + q);
    
    s[ 0 ] = u + v;
    num = 1;
  }
  
  // resubstitute 
  sub = 1.0/3 * A;
  
  for (i = 0; i < num; ++i)
    s[ i ] -= sub;
  
  return num;
}


G4int G4ToroidalSurface::SolveQuadric(G4double c[], G4double s[] )
{
  // From Graphics Gems I by Jochen Schwartz
  G4double p, q, D;
  
  // Normal form: x^2 + px + q = 0 
  p = c[ 1 ] / (2 * c[ 2 ]);
  q = c[ 0 ] / c[ 2 ];
  
  D = p * p - q;
  
  if (IsZero(D))
  {
    s[ 0 ] = - p;
    return 1;
  }
  else if (D < 0)
  {
    return 0;
  }
  else if (D > 0)
  {
    G4double sqrt_D = sqrt(D);
    
    s[ 0 ] =   sqrt_D - p;
    s[ 1 ] = - sqrt_D - p;
    return 2;
  }

  return 0;
}
