// 
// E.Medernach 2000
//

#define EPSILON 1e-12
#define INFINITY 1e+12
#define TORUSPRECISION 0.001 //1.0  // or whatever you want for precision (it is TorusEquation related)

#define NBPOINT 6
#define ITERATION 8 //20 But 8 is really enough for Newton with a good guess
#define NOINTERSECTION -1//kInfinity

#define DEBUGTORUS 0

/*
  Torus implementation with Newton Method and Bounding volume
 */


#define G4double double


#include <stdio.h>
#include <math.h>
#include "torus.h"

double cos(double x);
double sin(double x);


double sqrt(double x);
double fabs(double x);


inline int CheckAngle (double x,double y,double phi,double deltaphi)
{
  /** Note: this is possble to avoid atan by projecing -PI;PI to -inf;inf **/
  
		double theta ;
 
		theta = atan(x/y);
		if (y < 0.0) theta += M_PI;
		if (theta < 0.0) theta += 2*M_PI;

		if ((theta >= phi) && (theta <= (phi + deltaphi))) {
				return 1;
		} else {
				return 0;
		}
}

inline double IntersectPlanarSection (double x,double y,double dx,double dy,double phi,double deltaphi)
{
  /*** Intersect a ray with plan (phi) and (phi + deltaphi) ***/
  /*** the point is outside phi..phi+deltaphi ***/
  double Lambda1,Lambda2 ;
  Lambda1 = -(y - x*tan(phi))/(dy - dx*tan(phi));
  Lambda2 = -(y - x*tan(phi + deltaphi))/(dy - dx*tan(phi + deltaphi));
  if (Lambda1 < Lambda2) {
	return Lambda1;
  } else {
	return Lambda2;
  }
}

inline  double TorusEquation (x, y, z, R0, R1)
double x;
double y;
double z;
double R0;
double R1;
{
		/*
		  An interesting property is that the sign
		  tell if the point is inside or outside
		  or if > EPSILON on the surface
		*/
		double temp;

		temp = ((x*x + y*y + z*z) + R0*R0 - R1*R1) ;
		temp = temp*temp ;
		temp = temp - 4*R0*R0*(x*x + y*y) ;

		/*
		  > 0 Outside
		  < 0 Inside
		*/
		return temp ;
}


inline double TorusDerivativeX (x, y, z, R0, R1)
double x;
double y;
double z;
double R0;
double R1;
{
		return 4*x*(x*x + y*y + z*z +  R0*R0 - R1*R1) - 8*R0*R0*x ;
}

inline double TorusDerivativeY (x, y, z, R0, R1)
double x;
double y;
double z;
double R0;
double R1;
{
		return 4*y*(x*x + y*y + z*z +  R0*R0 - R1*R1) - 8*R0*R0*y ;
}


inline double TorusDerivativeZ (x, y, z, R0, R1)
double x;
double y;
double z;
double R0;
double R1;
{
		return 4*z*(x*x + y*y + z*z +  R0*R0 - R1*R1) ;
}


inline  double ParaboloidEquation (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return z - H*(x*x + y*y)/(L*L) ;
}

inline  double ParaboloidDerX (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return - 2*H*x/(L*L) ;
}

inline  double ParaboloidDerY (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return - 2*H*y/(L*L) ;
}

inline  double ParaboloidDerZ (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return 1 ;
}


inline  double HyperboloidEquation (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return (x*x + y*y) - z*z + H*H - L*L ;
}

inline  double HyperboloidDerX (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return 2*x ;
}

inline  double HyperboloidDerY (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return 2*y ;
}

inline  double HyperboloidDerZ (x, y, z, H, L)
     double x;
     double y;
     double z;
     double H;
     double L;
     
{
  return -2*z ;
}


void BVMParaboloidIntersection (G4double x,G4double y,G4double z,
				G4double dx,G4double dy,G4double dz,
				G4double H, G4double L,
				G4double *NewL,int *valid)
{
  /* We use the box [-L L]x[-L L]x[0 H] */
  /* there is only one interval at maximum */

  /* NewL and valid are array of 6 elements */

  if (dz != 0) {
    /* z = 0 */
    NewL[0] = -z/dz ;
    if ((fabs(x + NewL[0]*dx) < L) && (fabs(y + NewL[0]*dy) < L)) {
      valid[0] = 1;
    } else {
      valid[0] = 0;
    }
    
    /* z = H */
    NewL[1] = -(z-H)/dz ;
    if ((fabs(x + NewL[1]*dx) < L) && (fabs(y + NewL[1]*dy) < L)) {
      valid[1] = 1;
    } else {
      valid[1] = 0;
    }
    
  } else {
    NewL[0] = -1.0 ;
    NewL[1] = -1.0 ;
    valid[0] = 0;
    valid[1] = 0;
  }

  if (dx != 0) {
    /* x = -L */
    NewL[2] = -(x+L)/dx ;
    if ((fabs(z - H/2 +NewL[2]*dz) < H/2) && (fabs(y + NewL[2]*dy) < L)) {
      valid[2] = 1;
    } else {
      valid[2] = 0;
    }
    
    /* z = H */
    NewL[3] = -(x-L)/dx ;
    if ((fabs(z - H/2 + NewL[3]*dz) < H/2) && (fabs(y + NewL[3]*dy) < L)) {
      valid[3] = 1;
    } else {
      valid[3] = 0;
    }
    
  } else {
    NewL[2] = -1.0 ;
    NewL[3] = -1.0 ;
    valid[2] = 0;
    valid[3] = 0;
  }

  if (dy != 0) {
    /* y = -L */
    NewL[4] = -(y+L)/dy ;
    if ((fabs(z - H/2 +NewL[4]*dz) < H) && (fabs(y + NewL[4]*dy) < L)) {
      valid[4] = 1;
    } else {
      valid[4] = 0;
    }
    
    /* z = H */
    NewL[5] = -(y-L)/dy ;
    if ((fabs(z - H/2 + NewL[5]*dz) < H) && (fabs(y + NewL[5]*dy) < L)) {
      valid[5] = 1;
    } else {
      valid[5] = 0;
    }
    
  } else {
    NewL[4] = -1.0 ;
    NewL[5] = -1.0 ;
    valid[4] = 0;
    valid[5] = 0;
  }

}

void BVMHyperboloidIntersection (G4double x,G4double y,G4double z,
				G4double dx,G4double dy,G4double dz,
				G4double H, G4double L,
				G4double *NewL,int *valid)
{
  /* We use the box [-L L]x[-L L]x[-H H] */
  /* there is only one interval at maximum */

  /* NewL and valid are array of 6 elements */

  if (dz != 0) {
    /* z = -H */
    NewL[0] = -(z+H)/dz ;
    if ((fabs(x + NewL[0]*dx) < L) && (fabs(y + NewL[0]*dy) < L)) {
      valid[0] = 1;
    } else {
      valid[0] = 0;
    }
    
    /* z = H */
    NewL[1] = -(z-H)/dz ;
    if ((fabs(x + NewL[1]*dx) < L) && (fabs(y + NewL[1]*dy) < L)) {
      valid[1] = 1;
    } else {
      valid[1] = 0;
    }
    
  } else {
    NewL[0] = -1.0 ;
    NewL[1] = -1.0 ;
    valid[0] = 0;
    valid[1] = 0;
  }

  if (dx != 0) {
    /* x = -L */
    NewL[2] = -(x+L)/dx ;
    if ((fabs(z +NewL[2]*dz) < H) && (fabs(y + NewL[2]*dy) < L)) {
      valid[2] = 1;
    } else {
      valid[2] = 0;
    }
    
    /* z = H */
    NewL[3] = -(x-L)/dx ;
    if ((fabs(z + NewL[3]*dz) < H) && (fabs(y + NewL[3]*dy) < L)) {
      valid[3] = 1;
    } else {
      valid[3] = 0;
    }
    
  } else {
    NewL[2] = -1.0 ;
    NewL[3] = -1.0 ;
    valid[2] = 0;
    valid[3] = 0;
  }

  if (dy != 0) {
    /* y = -L */
    NewL[4] = -(y+L)/dy ;
    if ((fabs(z +NewL[4]*dz) < H) && (fabs(y + NewL[4]*dy) < L)) {
      valid[4] = 1;
    } else {
      valid[4] = 0;
    }
    
    /* z = H */
    NewL[5] = -(y-L)/dy ;
    if ((fabs(z + NewL[5]*dz) < H) && (fabs(y + NewL[5]*dy) < L)) {
      valid[5] = 1;
    } else {
      valid[5] = 0;
    }
    
  } else {
    NewL[4] = -1.0 ;
    NewL[5] = -1.0 ;
    valid[4] = 0;
    valid[5] = 0;
  }

}

void BVMIntersection(G4double x,G4double y,G4double z,
			      G4double dx,G4double dy,G4double dz,
			      G4double Rmax, G4double Rmin,
			      G4double *NewL,int *valid)
{

  if (dz != 0) {
    G4double DistToZ ;
    /* z = + Rmin */
    NewL[0] = (Rmin - z)/dz ;
    /* z = - Rmin */
    NewL[1] = (-Rmin - z)/dz ;
    /* Test validity here (*** To be optimized ***) */
    if (NewL[0] < 0.0) valid[0] = 0;
    if (NewL[1] < 0.0) valid[1] = 0;
    DistToZ = (x+NewL[0]*dx)*(x+NewL[0]*dx) + (y+NewL[0]*dy)*(y+NewL[0]*dy);
    if (DistToZ	- (Rmax + Rmin)*(Rmax + Rmin) > 0)
      valid[0] = 0;
    if (DistToZ	- (Rmax - Rmin)*(Rmax - Rmin) < 0)
      valid[0] = 0;
    DistToZ = (x+NewL[1]*dx)*(x+NewL[1]*dx) + (y+NewL[1]*dy)*(y+NewL[1]*dy);
    if (DistToZ	- (Rmax + Rmin)*(Rmax + Rmin) > 0)
      valid[1] = 0;
    if (DistToZ	- (Rmax - Rmin)*(Rmax - Rmin) < 0)
      valid[1] = 0;
  } else {
    /* if dz == 0 we could know the exact solution */
    /* Well, this is true but we have not expected precision issue from sqrt .. */
    NewL[0] = -1.0;
    NewL[1] = -1.0;
    valid[0] = 0;
    valid[1] = 0;
  }

  /* x� + y� = (Rmax + Rmin)� */
  if ((dx != 0) || (dy != 0)) {
    G4double a,b,c,d;
    
    a = dx*dx + dy*dy ;
    b = 2*(x*dx + y*dy) ;
    c = x*x + y*y - (Rmax + Rmin)*(Rmax + Rmin) ;
    d = b*b - 4*a*c ;
    
    if (d < 0) {
      valid[2] = 0;
      valid[3] = 0;
      NewL[2] = -1.0;
      NewL[3] = -1.0;
    } else {
      d = sqrt(d) ;
      NewL[2] = (d - b)/(2*a);
      NewL[3] = (-d - b)/(2*a);
      if (NewL[2] < 0.0) valid[2] = 0;
      if (fabs(z + NewL[2]*dz) - Rmin > EPSILON) valid[2] = 0;
      if (NewL[3] < 0.0) valid[3] = 0;
      if (fabs(z + NewL[3]*dz) - Rmin > EPSILON) valid[3] = 0;
    }
  } else {
    /* only dz != 0 so we could know the exact solution */
    /* this depends only for the distance to Z axis */
    /* BUT big precision problem near the border.. */
    /* I like so much Newton to increase precision you know.. => */

    NewL[2] = -1.0;
    NewL[3] = -1.0;
    valid[2] = 0;
    valid[3] = 0;
	
    /*** Try This to see precision issue with sqrt(~ 0)
	 G4double DistToZ ;
	 G4double result;
	 G4double guess;
	
	 DistToZ = sqrt(x*x + y*y) ;
	
	 if ((DistToZ < (Rmax - Rmin)) || (DistToZ > (Rmax + Rmin))) {
	 return -1.0 ;
	 }
	
	 result = sqrt((Rmin + Rmax - DistToZ)*(Rmin - Rmax + DistToZ));

	 if (dz < 0) {
	 if (z > result) {
	 return (result - z)/dz;
	 } else {
	 if (z > -result) {
	 return (-result - z)/dz;
	 } else 
	 return -1.0;
	 }
	 } else {
	 if (z < -result) {
	 return (z + result)/dz;
	 } else {
	 if (z < result) {
	 return (z - result)/dz;
	 } else 
	 return -1.0;
	 }
	 }
    */
  }
  

  /* x� + y� = (Rmax - Rmin)� */
  if ((dx != 0) || (dy != 0)) {
    G4double a,b,c,d;
    
    a = dx*dx + dy*dy ;
    b = 2*(x*dx + y*dy) ;
    c = x*x + y*y - (Rmax - Rmin)*(Rmax - Rmin) ;
    d = b*b - 4*a*c ;
    
    if (d < 0) {
      valid[4] = 0;
      valid[5] = 0;
      NewL[4] = -1.0;
      NewL[5] = -1.0;
    } else {
      d = sqrt(d) ;
      NewL[4] = (d - b)/(2*a);
      NewL[5] = (-d - b)/(2*a);
      if (NewL[4] < 0.0) valid[4] = 0;
      if (fabs(z + NewL[4]*dz) - Rmin > EPSILON) valid[4] = 0;
      if (NewL[5] < 0.0) valid[5] = 0;
      if (fabs(z + NewL[5]*dz) - Rmin > EPSILON) valid[5] = 0;
    }
  } else {
    /* only dz != 0 so we could know the exact solution */
    /* OK but same as above .. */
    valid[4] = 0;
    valid[5] = 0;
    NewL[4] = -1.0;
    NewL[5] = -1.0;
  }
}

void SortIntervals (int NbElem,G4double *SortL,G4double *NewL,int *valid,int *NbIntersection)
{
  int i,j;
  G4double swap;
	
  (*NbIntersection) = 0;
  SortL[0] = -INFINITY;
	
  for (i=0;i<NbElem;i++) {
    if (valid[i] != 0) {
      SortL[(*NbIntersection)] = NewL[i] ;
      for (j=(*NbIntersection);j>0;j--) {
	if (SortL[j] < SortL[j-1]) {
	  swap = SortL[j-1] ;
	  SortL[j-1] = SortL[j];
	  SortL[j] = swap;
	}
      }
		
      (*NbIntersection) ++;
    }
  }
  /* Delete double value */
  /* When the ray hits a corner we have a double value */
  for (i=0;i<(*NbIntersection)-1;i++) {
    if (SortL[i+1] - SortL[i] < EPSILON) {
      if (((*NbIntersection) & (1)) == 1) {
	/* If the NbIntersection is odd then we keep one value */
	for (j=i+1;j<(*NbIntersection);j++) {
	  SortL[j-1] = SortL[j] ;
	}
	(*NbIntersection) --;
      } else {
	/* If it is even we delete the 2 values */
	for (j=i+2;j<(*NbIntersection);j++) {
	  SortL[j-2] = SortL[j] ;
	}
	(*NbIntersection) -= 2;
      }
    }
  }
}


/* TODO:
  check if the root is entering the torus (with gradient)
  clean problems when Rmin ~ Rmax (BVM is not good when near Z axis)
 */

/** Now the interesting part .. **/

int SafeNewton(G4double x, G4double y, G4double z,
		    G4double dx, G4double dy, G4double dz,
		    G4double Rmax, G4double Rmin,
		    G4double *Lmin,G4double *Lmax)
{
  /** SafeNewton is a clipping interval Newton method **/
  G4double P[5][2],D[2] ;
  G4double Lx,Ly,Lz ;
  G4double NewMin,NewMax;
  
  int IntervalIsVoid = 1;
  int NewtonIsSafe = 0;
  
  /*** Calculating Control Points  ***/
  
  /*
    0    	p0 = F((*Lmin))
    1/4  	p1 = F((*Lmin)) + ((*Lmax) - (*Lmin))/4 * F'((*Lmin))
    2/4  	p2 = 1/6 * (32*F(((*Lmax) + (*Lmin))/2) - (p0 + 4*p1 + 4*p3 + p4))  
    3/4         p3 = F((*Lmax)) - ((*Lmax) - (*Lmin))/4 * F'((*Lmax))
    1           p4 = F((*Lmax))
  */

  
  Lx = x + (*Lmin)*dx;
  Ly = y + (*Lmin)*dy;
  Lz = z + (*Lmin)*dz;

  D[0] = dx*HyperboloidDerX(Lx,Ly,Lz,Rmax,Rmin);
  D[0] += dy*HyperboloidDerY(Lx,Ly,Lz,Rmax,Rmin);
  D[0] += dz*HyperboloidDerZ(Lx,Ly,Lz,Rmax,Rmin);

  P[0][0] = (*Lmin);
  P[0][1] = HyperboloidEquation(Lx,Ly,Lz,Rmax,Rmin);

  if (fabs(P[0][1]) < TORUSPRECISION) {
    NewtonIsSafe = 1;
    //fprintf(stderr,"(fabs(P[0][1]) < TORUSPRECISION)\n");
    return NewtonIsSafe;
  }
  
  if (((*Lmax) - (*Lmin)) < EPSILON) {
    //fprintf(stderr,"(((*Lmax) - (*Lmin)) < EPSILON)\n");
    return 1;
  }

  P[1][0] = (*Lmin) + ((*Lmax) - (*Lmin))/4;
  P[1][1] = P[0][1] + (((*Lmax) - (*Lmin))/4.0) * D[0];
  
  Lx = x + (*Lmax)*dx;
  Ly = y + (*Lmax)*dy;
  Lz = z + (*Lmax)*dz;

  D[1] = dx*HyperboloidDerX(Lx,Ly,Lz,Rmax,Rmin);
  D[1] += dy*HyperboloidDerY(Lx,Ly,Lz,Rmax,Rmin);
  D[1] += dz*HyperboloidDerZ(Lx,Ly,Lz,Rmax,Rmin);

  P[4][0] = (*Lmax);
  P[4][1] = HyperboloidEquation(Lx,Ly,Lz,Rmax,Rmin);
  P[3][0] = (*Lmax) - ((*Lmax) - (*Lmin))/4;
  P[3][1] = P[4][1] - ((*Lmax) - (*Lmin))/4 * D[1];

  Lx = x + ((*Lmax)+(*Lmin))/2*dx;
  Ly = y + ((*Lmax)+(*Lmin))/2*dy;
  Lz = z + ((*Lmax)+(*Lmin))/2*dz;

  P[2][0] = ((*Lmax) + (*Lmin))/2;
  P[2][1] = (16*HyperboloidEquation(Lx,Ly,Lz,Rmax,Rmin) - (P[0][1] + 4*P[1][1] + 4*P[3][1] + P[4][1]))/6 ;

  
  
  //fprintf(stderr,"\n");
  //fprintf(stderr,"Lmin = %14f\n",(*Lmin));
  //fprintf(stderr,"Lmax = %14f\n",(*Lmax));
  //fprintf(stderr,"P[0] = %14f\n",P[0][1]);
  //fprintf(stderr,"P[1] = %14f\n",P[1][1]);
  //fprintf(stderr,"P[2] = %14f\n",P[2][1]);
  //fprintf(stderr,"P[3] = %14f\n",P[3][1]);
  //fprintf(stderr,"P[4] = %14f\n",P[4][1]);

#if DEBUGTORUS
  G4cout << "G4Torus::SafeNewton    Lmin = " << (*Lmin) << G4endl ;
  G4cout << "G4Torus::SafeNewton    Lmax = " << (*Lmax) << G4endl ;
  G4cout << "G4Torus::SafeNewton    P[0] = " << P[0][1] << G4endl ;
  G4cout << "G4Torus::SafeNewton    P[1] = " << P[1][1] << G4endl ;
  G4cout << "G4Torus::SafeNewton    P[2] = " << P[2][1] << G4endl ;
  G4cout << "G4Torus::SafeNewton    P[3] = " << P[3][1] << G4endl ;
  G4cout << "G4Torus::SafeNewton    P[4] = " << P[4][1] << G4endl ;
#endif

  /** Ok now we have all control points, we could compute the convex area **/
  /** Problems:
      - if there is one point with a ~ 0 coordinate and all the other the same sign we
      miss the value
      - if there are more than a root in the interval then the interval length does not
      decrease to 0. A solution may be to split intervals in the middle but how to
      know that we must split ?
      - we have to compute convex area of the control point before applying intersection
      with y=0
  **/

  /*** For each points make 2 sets. A set of positive points and a set of negative points ***/
  /*** Note: could be better done with scalar product .. ***/

  /* there is an intersection only if each have different signs */
  /* PROBLEM : If a control point have a 0.00 value the sign check is wrong */
  {
    G4double Intersection ;
    int i,j;

    NewMin = (*Lmax) ;
    NewMax = (*Lmin) ;

    for (i=0;i<5;i++)
      for (j=i+1;j<5;j++)
	{
	  /* there is an intersection only if each have different signs */
	  if (((P[j][1] > -TORUSPRECISION) && (P[i][1] < TORUSPRECISION)) ||
	      ((P[j][1] < TORUSPRECISION) && (P[i][1] > -TORUSPRECISION))) {
	    IntervalIsVoid  = 0;
	    Intersection = P[j][0] - P[j][1]*((P[i][0] - P[j][0])/(P[i][1] - P[j][1]));
	    if (Intersection < NewMin) {
	      NewMin = Intersection;
	    }
	    if (Intersection > NewMax) {
	      NewMax = Intersection;
	    }
	  }
	}
    if (IntervalIsVoid != 1) {
      (*Lmax) = NewMax;
      (*Lmin) = NewMin;
    }
  }
  
  if (IntervalIsVoid == 1) {
    //fprintf(stderr,"(IntervalIsVoid == 1)\n");
    return -1;
  }
  
  //fprintf(stderr,"NewMin = %f  NewMax = %f\n",NewMin,NewMax);
  /** Now we have each Extrema point of the new interval **/
  
  return NewtonIsSafe;
}


G4double Newton (G4double guess,
		 G4double x, G4double y, G4double z,
		 G4double dx, G4double dy, G4double dz,
		 G4double Rmax, G4double Rmin,
		 G4double Lmin,G4double Lmax)
{
  /* So now we have a good guess and an interval where if there are an intersection the root must be */

  G4double Lx = 0;
  G4double Ly = 0;
  G4double Lz = 0;
  G4double Value = 0;
  G4double Gradient = 0;
  G4double Lambda ;

  int i=0;

  /* Reduce interval before applying Newton Method */
  {
    int NewtonIsSafe ;

    while ((NewtonIsSafe = SafeNewton(x,y,z,dx,dy,dz,Rmax,Rmin,&Lmin,&Lmax)) == 0) ;

    guess = Lmin;
  }

  /*** BEWARE ***/
  /* A typical problem is when Gradient is zero */
  /* This is due to some 0 values in point or direction */
  /* To solve that we move a little the guess 
  if ((((x == 0) || (y == 0)) || (z == 0)) || 
      (((dx == 0) || (dy == 0)) || (dz == 0))) 
      guess += EPSILON;*/
  
  Lambda = guess;
  Value = HyperboloidEquation(x + Lambda*dx,y + Lambda*dy,z + Lambda*dz,Rmax,Rmin);

  //fprintf(stderr,"NEWTON begin with L = %f and V = %f\n",Lambda,Value);
  
  /*** Beware: we must eliminate case with no root ***/
  /*** Beware: In some rare case we converge to the false root (internal border)***/
  /***
   {
    FILE *fi;
    int i;
    fi = fopen("GNUplot.out","w+");
    //fprintf(fi,"# Newton plot\n");
	  
    for (i = 0; i < 1000 ; i ++) {
      Lx = x + (Lmin + i*(Lmax - Lmin)/1000.0)*dx;
      Ly = y + (Lmin + i*(Lmax - Lmin)/1000.0)*dy;
      Lz = z + (Lmin + i*(Lmax - Lmin)/1000.0)*dz;
      Value = HyperboloidEquation(Lx,Ly,Lz,Rmax,Rmin);
      //fprintf(fi," %f %f\n",Lmin + i*(Lmax - Lmin)/1000.0,Value );
    }
	  
    fclose(fi);
  }
   
  ***/
      
  /* In fact The Torus Equation give big number so TORUS PRECISION is not EPSILON */	  
  while (/* ?? (fabs(Value/Gradient) > 1e-2) ||*/ (fabs(Value) > TORUSPRECISION)) {
	  
    //  do {
    Lx = x + Lambda*dx;
    Ly = y + Lambda*dy;
    Lz = z + Lambda*dz;
    Value = HyperboloidEquation(Lx,Ly,Lz,Rmax,Rmin);

    Gradient = dx*HyperboloidDerX(Lx,Ly,Lz,Rmax,Rmin);
    Gradient += dy*HyperboloidDerY(Lx,Ly,Lz,Rmax,Rmin);
    Gradient += dz*HyperboloidDerZ(Lx,Ly,Lz,Rmax,Rmin);

	/*
	if (Gradient > -EPSILON)
	  return Lmin;
	*/
	
    /***
	if ((beware != 0) && (Gradient > -EPSILON)) {
    ***/

    /** Newton does not go to the root because interval is too big **/
    /** In fact Newton is known to converge if |f.f''/(f'^2)| < 1 **/
    /** There is two cases: ray hits or not **/
    /** If ray hits we must search for a better intervals **/
    /** but if there are no hits then we could not .. **/
    /** So the easier way the best: if Newton encounter a problem
	it says to the BVM that the guess is no good
	then the BVM search for a better intervals, possibly none
	in this case no intersection, else we go back to Newton 
    **/

    /**
       Perhaps we have not to break Newton at the beginning because we could converge after some move
       May be not: If we are here this means that the root we want is rejecting. We could converge to 
       another root.
       PROBLEMS
    **/
    /* root is repulsive from this guess could you give me another guess ? 
       Note: that it may be no root in this area ..
       Note: Lmin and Lmax are always outside the torus as a part of the BVM.
       We just want a point in this direction with a gradient < 0
				 
       guess = FindABetterGuess(Rmax,Rmin,guess,Lmin,Lmax);
    */
    Lambda = Lambda - Value/Gradient ;

#if DEBUGTORUS
    G4cout << "Newton     Iteration " << i << G4endl ;
    G4cout << "Newton     Lambda = " << Lambda << " Value = " << Value << " Grad = " << Gradient << G4endl;
    G4cout << "Newton     Lmin = " << Lmin << " Lmax = " << Lmax << G4endl ;
#endif
    //fprintf(stderr,"Newton    Iteration %d\n",i);
    //fprintf(stderr,"Newton    Lambda = %f  Value = %f  Grad = %f\n",Lambda,Value,Gradient);
    
    i ++;

    if (i > ITERATION) 
      return NOINTERSECTION; //no convergency ??

  } //while (/* ?? (fabs(Value/Gradient) > 1e-2) ||*/ (fabs(Value) > TORUSPRECISION));


#if DEBUGTORUS
  G4cout << "Newton    Exiting with Lambda = " << Lambda << G4endl ;
  G4cout << "Newton    Exiting with Value = " << Value << G4endl ;
#endif

  //just a check 
  if (Lambda < 0.0) {
	//fprintf(stderr,"Newton end with a negative solution ..\n");
	return NOINTERSECTION;
  }
  //fprintf(stderr,"NEWTON: Lamdba = %f\n",Lambda);
  return Lambda ;
}

/*
G4double DistanceToTorus (G4double x,G4double y,G4double z,
				   G4double dx,G4double dy,G4double dz,
				   G4double Rmax,G4double Rmin)
*/
double DistanceToTorus (Intersect * Inter)
{
  static int Vstatic = 0;
  G4double Lmin,Lmax;
  G4double guess;
  G4double SortL[4];
   
  int NbIntersection = 0;

  G4double NewL[NBPOINT];
  int valid[] = {1,1,1,1,1,1} ;
  int j;

  		double x,y,z,dx,dy,dz;
		double Rmax,Rmin;
		double phi,deltaphi;

		j = 0;


  		dx = Inter->dx;
		dy = Inter->dy;
		dz = Inter->dz;
		x = Inter->x;
		y = Inter->y;
		z = Inter->z;
		Rmax = Inter->R0 ;
		Rmin = Inter->R1 ;
		phi = Inter->phi;
		deltaphi = Inter->deltaphi;


  /*** Compute Intervals  from Bounding Volume ***/

		//BVMIntersection(x,y,z,dx,dy,dz,Rmax,Rmin,NewL,valid);
		BVMHyperboloidIntersection(x,y,z,dx,dy,dz,Rmax,Rmin,NewL,valid);

  /*
    We could compute intervals value 
    Sort all valid NewL to SortL.
    There must be 4 values at max and 
    odd one if point is inside
  */

  SortIntervals(6,SortL,NewL,valid,&NbIntersection);
  if (BVM_ONLY == 1)
    return SortL[0] ; 

#if 0
  // Torus Only
  {
    /*** Length check ***/
    G4double LengthMin = 0.82842712*Rmin;
				
    switch(NbIntersection) {
    case 1:
      if (SortL[0] < EPSILON) {
	if (fabs(HyperboloidEquation(x,y,z,Rmax,Rmin)) < TORUSPRECISION) {
	  return 0.0;
	} else {
	  return NOINTERSECTION;
	}
      }
      break;
    case 2:
      if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection = 0;
      break;
    case 3:
      if (SortL[0] < EPSILON) {
	if (fabs(HyperboloidEquation(x,y,z,Rmax,Rmin)) < TORUSPRECISION) {
	  return 0.0;
	} else {
	  NbIntersection --;
	  SortL[0] = SortL[1] ;
	  SortL[1] = SortL[2] ;
	  if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection = 0;
	}
      } else {
	if ((SortL[2] - SortL[1]) < LengthMin) NbIntersection -= 2;
      }
      break;
    case 4:
      if ((SortL[1] - SortL[0]) < LengthMin) {
	NbIntersection -= 2;
	SortL[0] = SortL[2];
	SortL[1] = SortL[3];
	if ((SortL[1] - SortL[0]) < LengthMin) NbIntersection -= 2;	
      }
      break;
    }
  }
#endif
  
#if DEBUGTORUS
  {
    int i;
    G4cout.precision(16);
    G4cout << "DistanceToTorus    INTERVALS" << G4endl ;
    for (i=0;i<NbIntersection;i++) {
      G4cout << "DistanceToTorus    " << SortL[i] << G4endl ;
    }
  }
#endif

  Vstatic ++;
   
  //if  ((Vstatic % 2) == 0) return SortL[0];
  //printf("NbIntersection = %d\n",NbIntersection);
  
  
  /* BVM Test 

     switch(NbIntersection) {
     case 0:
     return -1.0;
     break;
     case 1:
     return -1.0;
     break;
     case 2:
     return -1.0;
     break;
     case 3:
     return -1.0;
     break;
     case 4:
     return -1.0;
     break;
     }
  */
  			
  /*** If the ray intersects the torus it necessary intersects the BVMax ***/
  /*** So it is necessary into *an* interval from the BVM ***/

  /** Note : In general there are only 2 intersections so computing the second interval
      could be done only if the first one does not contain any root */

  /* NOW there is 2 possibilities */
  /* If inside the BVM (or Torus instead), take "0, SortL[0] .." */
  /* If outside the BVM, we have intervals where if there is an intersection the root must be */
  /* Now Lmin1 <= Lambda <= Lmax and there is a *unique* root */
  /* Newton Methods in this interval from the guess */

  /*** Beware The first interval could be the bad one and we have to see other one ***/
  /*** We must have a way to decide if an interval contains root or not .. ***/

  /*** 
       Beware: If the original point is near the torus (into the BVM not the torus)
       we have serious precision issue (bad guess value) try it with a big Rmin
  ***/

  /* We are Inside the BVM if the number of intersection is odd */
  /* Not necessary an intersection with Torus if point outside Torus and Inside BVM ! */

  if (((NbIntersection) & (1)) != 0) {
    /*** If we are Inside the BVM Lmin = 0. Lmax is the point ***/
    /***    there is necessary an intersection if the point is inside the Torus ***/
    int InsideTorus = 0;

    Lmin = 0.0 ;
    Lmax  = SortL[0] ;

    if (HyperboloidEquation(x,y,z,Rmax,Rmin) < 0.0) {
      
      InsideTorus = 1;
      /* As we are inside the torus it must have an intersection */
      /* To have a good guess we take Lmax - Rmin/8.0 */
      /*(What is the best value for a square to be like a circle ?) */
      /* If we are inside the torus the upper bound is better */
	  //return 1000.0;
      guess =  Lmax - Rmin*0.125;
      //printf("DistanceToTorus    Inside the torus\n");
      
#if DEBUGTORUS
      G4cout << "DistanceToTorus    Inside the torus" << G4endl ;
      G4cout << "DistanceToTorus    Initial Guess is " << guess << G4endl ;
#endif
      
    } else {
	  //	  return 1000.0;
      //printf("DistanceToTorus    Outside the torus\n");
#if DEBUGTORUS
      G4cout.precision(16);
      G4cout << "DistanceToTorus    point " << x << ", " << y << ", " << z << ", "  << " is outside the torus "
	     << " Rmax = " << Rmax << " Rmin = " << Rmin << " Teq = " << HyperboloidEquation(x,y,z,Rmax,Rmin) << G4endl ;
#endif
      InsideTorus = 0;
      /* PROBLEMS what to choose ? 0.0 ? */
      /* 0.0 is generally a good guess, but there is case that it is very bad (hit center torus when inside BVM) */

      if (Lmax > Rmin) {
	/* we are in the case where we hit center torus */
	
	//return 100000.0;
	guess = Lmax;
	
      } else {
	/* general case */
	guess = 0.0;
      }
    }

    /* Ready to do Newton */
    guess = Newton(guess,x,y,z,dx,dy,dz,Rmax,Rmin,Lmin,Lmax);

#if DEBUGTORUS
    G4cout << "DistanceToTorus    First Newton guess = " << guess << G4endl ;
    G4cout << "DistanceToTorus    Lmin = " << Lmin << "  Lmax = " << Lmax << G4endl ;
#endif

    /* In case we are the origin point is just in the surface 
       the NbIntersection will be odd and guess will be zero
       Anyway, it is correct to say that distance is zero but
       we want to return +inf if we are exiting the solid
       So ..
    */

    /* Check here is the root found is into interval */

    if ((guess >= (Lmin - EPSILON)) && (guess <= (Lmax + EPSILON))) {
      return guess ;
    } else {
      if (NbIntersection == 3) {
				/** OK we are in the small part around the BVM **/
				/** So we check the second interval **/
	Lmin = SortL[1];
	Lmax = SortL[2];
	guess = Lmin;
	
	guess = Newton(guess,x,y,z,dx,dy,dz,Rmax,Rmin,Lmin,Lmax);
#if DEBUGTORUS
	G4cout << "DistanceToTorus    Second Newton guess = " << guess << G4endl ;
	G4cout << "DistanceToTorus    Lmin = " << Lmin << "  Lmax = " << Lmax << G4endl ;
#endif
	if ((guess >= (Lmin - EPSILON)) && (guess <= (Lmax + EPSILON))) {
	  return guess;
	} else {
	  return NOINTERSECTION;
	}
      } else {
	if (InsideTorus == 1) {
	  /* Incredible : sometimes precisions errors bring us here 
	     with guess = SortL[0]
	     So we return guess ..

	     PROBLEMS 99%
	  
										
	  printf("Torus: Root not found final (guess - Limit) = %f\n"
		 ,guess - SortL[0]);
	  printf("point: %f %f %f\n",x,y,z);
	  printf("dir  : %f %f %f\n",dx,dy,dz);
	  */
	  
	  return 100000.0;//guess;
	  exit(1);
										
	}
	return NOINTERSECTION;
      }
    }
	
	

  } else { // Outside
    /*** If we are Out then we need more to know if intersection exists ***/
    /***  there is 2 intersection points at least (perhaps the same) with BVMax ***/

    /*** Return if no intersection with BVMax ***/

    if (NbIntersection == 0) 
      return NOINTERSECTION ;
				

    Lmin = SortL[0] ;
    Lmax = SortL[1] ;
    /** Lmin because it is probably near the BVM entry point **/
    /** PROBLEM if the ray hits the top of BVM with a small angle
	then the interval is too big and the guess is bad **/

    guess = Lmin ; 
				

    /*** We know only that if there is a solution, it is between Lmin and Lmax ***/
    /*** But we are not sure that there is one ... ***/
	
    /* Ready to do Newton */
    guess = Newton(guess,x,y,z,dx,dy,dz,Rmax,Rmin,Lmin,Lmax);

#if DEBUGTORUS
    G4cout << "DistanceToTorus    Newton with 2 or 4 points : " << guess << G4endl ;
#endif    

    /* Check here is the root found is into interval */
    if ((guess >= (Lmin - EPSILON)) && (guess <= (Lmax + EPSILON))) {
#if DEBUGTORUS
    G4cout << "DistanceToTorus    Newton gives a point into interval (Ok)" << G4endl ;
#endif    
      return guess;
    } else { 				
#if DEBUGTORUS
    G4cout << "DistanceToTorus    Newton does not give a point into interval (Ko)" << G4endl ;
#endif    
      if (NbIntersection == 4) {
	/* Well if that does not converge with the first interval try with the other one */
	Lmin = SortL[2] ;
	Lmax = SortL[3] ;

	guess = Lmin;
	guess = Newton(guess,x,y,z,dx,dy,dz,Rmax,Rmin,Lmin,Lmax);
	
	if ((guess >= (Lmin - EPSILON)) && (guess <= (Lmax + EPSILON))) {
	  return guess;
	} else {
	  return NOINTERSECTION;
	}
      } else {
	/* Certainly this is due to the BVM part that is not in Torus */

	return NOINTERSECTION ;
      }
    }
  }
}

inline G4double TorusGradient(G4double dx,
				       G4double dy,
				       G4double dz,
				       G4double x,
				       G4double y,
				       G4double z,
				       G4double Rmax,
				       G4double Rmin)
{
  /* This tell the normal at a surface point */
  G4double result;
  result = 0;
  result += dx*HyperboloidDerX(x,y,z,Rmax,Rmin); 
  result += dy*HyperboloidDerY(x,y,z,Rmax,Rmin); 
  result += dz*HyperboloidDerZ(x,y,z,Rmax,Rmin); 

  return result;
}


inline G4double ParaboloidGradient(G4double dx,
				       G4double dy,
				       G4double dz,
				       G4double x,
				       G4double y,
				       G4double z,
				       G4double Rmax,
				       G4double Rmin)
{
  /* This tell the normal at a surface point */
  G4double result;
  result = 0;
  result += dx*ParaboloidDerX(x,y,z,Rmax,Rmin); 
  result += dy*ParaboloidDerY(x,y,z,Rmax,Rmin); 
  result += dz*ParaboloidDerZ(x,y,z,Rmax,Rmin); 

  return result;
}

inline G4double HyperboloidGradient(G4double dx,
				       G4double dy,
				       G4double dz,
				       G4double x,
				       G4double y,
				       G4double z,
				       G4double Rmax,
				       G4double Rmin)
{
  /* This tell the normal at a surface point */
  G4double result;
  result = 0;
  result += dx*HyperboloidDerX(x,y,z,Rmax,Rmin); 
  result += dy*HyperboloidDerY(x,y,z,Rmax,Rmin); 
  result += dz*HyperboloidDerZ(x,y,z,Rmax,Rmin); 

  return result;
}

