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
// G4ToroidalSurface.cc
//
// ----------------------------------------------------------------------

#include "G4ToroidalSurface.hh"
#include "G4PhysicalConstants.hh"

G4ToroidalSurface::G4ToroidalSurface()
 : MinRadius(0.), MaxRadius(0.), TransMatrix(0), EQN_EPS(1e-9)
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
  delete TransMatrix;
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


G4Vector3D G4ToroidalSurface::SurfaceNormal(const G4Point3D&) const 
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
  G4Point3D  Base = Ray.GetStart();   // Base of the intersection ray
  G4Vector3D DCos = Ray.GetDir();     // Direction cosines of the ray
  G4int	     nhits=0;		      // Number of intersections
  G4double   rhits[4];		      // Intersection distances
  G4double   hits[4] = {0.,0.,0.,0.}; // Ordered intersection distances
  G4double   rho, a0, b0;	      // Related constants		
  G4double   f, l, t, g1, q, m1, u;   // Ray dependent terms		
  G4double   C[5];		      // Quartic coefficients	      
	
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
  g1 = f + rho * DCos.y()*DCos.y();
  q = a0 / (g1*g1);
  m1 = (l + 2.*rho*DCos.y()*Base.y()) / g1;
  u = (t +    rho*Base.y()*Base.y() + b0) / g1;
	
  //	Compute the coefficients of the quartic.			
  
  C[4] = 1.0;
  C[3] = 2. * m1;
  C[2] = m1*m1 + 2.*u - q*f;
  C[1] = 2.*m1*u - q*l;
  C[0] = u*u - q*t;
	
  //	Use quartic root solver found in "Graphics Gems" by Jochen Schwarze.
  nhits = SolveQuartic (C,rhits);
  
  //	SolveQuartic returns root pairs in reversed order.		
  m1 = rhits[0]; u = rhits[1]; rhits[0] = u; rhits[1] = m1;
  m1 = rhits[2]; u = rhits[3]; rhits[2] = u; rhits[3] = m1;

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
    // G4int Closest=0;
    
    //	    if(Counter>1)
    //		if(rhits[InsideBox[0]] > rhits[InsideBox[1]])
    //		    Closest=1;
    
    // Project polygon and do point in polygon
    // Projection also for curves etc.
    // Should probably be implemented in the curve class. 
    // G4Plane Plane1 = Ray.GetPlane(1);
    // G4Plane Plane2 = Ray.GetPlane(2);
    
    // Point in polygon
    return 1;
  }
  return 0;
}


G4int G4ToroidalSurface::SolveQuartic(G4double cc[], G4double ss[]  )
{
  // From Graphics Gems I by Jochen Schwartz
  
  G4double  coeffs[ 4 ];
  G4double  z, u, v, sub;
  G4double  A, B, C, D;
  G4double  sq_A, p, q, r;
  G4int     i, num;
  
    // Normal form: x^4 + Ax^3 + Bx^2 + Cx + D = 0 

  A = cc[ 3 ] / cc[ 4 ];
  B = cc[ 2 ] / cc[ 4 ];
  C = cc[ 1 ] / cc[ 4 ];
  D = cc[ 0 ] / cc[ 4 ];

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
    
    num = SolveCubic(coeffs, ss);
    
    ss[ num++ ] = 0;
  }
  else
  {
    // solve the resolvent cubic ... 
    coeffs[ 0 ] = 1.0/2 * r * p - 1.0/8 * q * q;
    coeffs[ 1 ] = - r;
    coeffs[ 2 ] = - 1.0/2 * p;
    coeffs[ 3 ] = 1;
    
    (void) SolveCubic(coeffs, ss);
    
    // ... and take the one real solution ... 
    z = ss[ 0 ];

    // ... to Build two quadric equations 
    u = z * z - r;
    v = 2 * z - p;
    
    if (IsZero(u))
      u = 0;
    else if (u > 0)
      u = std::sqrt(u);
    else
      return 0;

    if (IsZero(v))
      v = 0;
    else if (v > 0)
      v = std::sqrt(v);
    else
      return 0;

    coeffs[ 0 ] = z - u;
    coeffs[ 1 ] = q < 0 ? -v : v;
    coeffs[ 2 ] = 1;
    
    num = SolveQuadric(coeffs, ss);
    
    coeffs[ 0 ]= z + u;
    coeffs[ 1 ] = q < 0 ? v : -v;
    coeffs[ 2 ] = 1;
    
    num += SolveQuadric(coeffs, ss + num);
  }
  
  // resubstitute 
  
  sub = 1.0/4 * A;
  
  for (i = 0; i < num; ++i)
    ss[ i ] -= sub;
  
  return num;
}


G4int G4ToroidalSurface::SolveCubic(G4double cc[], G4double ss[]  )
{
  // From Graphics Gems I bu Jochen Schwartz
  G4int     i, num;
  G4double  sub;
  G4double  A, B, C;
  G4double  sq_A, p, q;
  G4double  cb_p, D;

  // Normal form: x^3 + Ax^2 + Bx + C = 0 
  A = cc[ 2 ] / cc[ 3 ];
  B = cc[ 1 ] / cc[ 3 ];
  C = cc[ 0 ] / cc[ 3 ];
  
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
      ss[ 0 ] = 0;
      num = 1;
    }
    else // one single and one G4double solution 
    {
      G4double u = std::pow(-q,1./3.);
      ss[ 0 ] = 2 * u;
      ss[ 1 ] = - u;
      num = 2;
    }
  }
  else if (D < 0) // Casus irreducibilis: three real solutions
  {
    G4double phi = 1.0/3 * std::acos(-q / std::sqrt(-cb_p));
    G4double t = 2 * std::sqrt(-p);
    
    ss[ 0 ] =   t * std::cos(phi);
    ss[ 1 ] = - t * std::cos(phi + pi / 3);
    ss[ 2 ] = - t * std::cos(phi - pi / 3);
    num = 3;
  }
  else // one real solution 
  {
    G4double sqrt_D = std::sqrt(D);
    G4double u = std::pow(sqrt_D - q,1./3.);
    G4double v = - std::pow(sqrt_D + q,1./3.);
    
    ss[ 0 ] = u + v;
    num = 1;
  }
  
  // resubstitute 
  sub = 1.0/3 * A;
  
  for (i = 0; i < num; ++i)
    ss[ i ] -= sub;
  
  return num;
}


G4int G4ToroidalSurface::SolveQuadric(G4double cc[], G4double ss[] )
{
  // From Graphics Gems I by Jochen Schwartz
  G4double p, q, D;
  
  // Normal form: x^2 + px + q = 0 
  p = cc[ 1 ] / (2 * cc[ 2 ]);
  q = cc[ 0 ] / cc[ 2 ];
  
  D = p * p - q;
  
  if (IsZero(D))
  {
    ss[ 0 ] = - p;
    return 1;
  }
  else if (D < 0)
  {
    return 0;
  }
  else if (D > 0)
  {
    G4double sqrt_D = std::sqrt(D);
    
    ss[ 0 ] =   sqrt_D - p;
    ss[ 1 ] = - sqrt_D - p;
    return 2;
  }

  return 0;
}
