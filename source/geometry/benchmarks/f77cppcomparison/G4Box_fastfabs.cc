// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Box_fastfabs.cc,v 1.1 1999-01-08 16:31:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Implementation for G4Box class
//

#include "G4Box.hh"

//#include "G4VoxelLimits.hh"
//#include "G4Transform.hh"

//#ifdef G4VISUALIZE
//#include "G4VWindow.hh"
//#include "G4Polyline.hh"
//#endif

//#include <math.h>

inline G4double fastfabs(const G4double p)
{
    return (p>=0) ? p : -p ;
}

// Private (implementation) enum: Not for external use   	
// Codes for faces (kPX=plus x face,kMY= minus y face etc)
enum ESide {kPX,kMX,kPY,kMY,kPZ,kMZ};

// Constructor - check & set half widths
G4Box::G4Box(const G4double pX,
             const G4double pY,
             const G4double pZ)
{
    fDx=pX; fDy=pY; fDz=pZ;
}


// Return whether point inside/outside/on surface, using tolerance
EInside G4Box::Inside(const G4ThreeVector& p) const
{
    EInside in=kOutside;

    if (fastfabs(p.x())<=fDx-kCarTolerance*0.5)
	{
	    if (fastfabs(p.y())<=fDy-kCarTolerance*0.5)
		{
		    if (fastfabs(p.z())<=fDz-kCarTolerance*0.5)
			{
			    in=kInside;
			}
		    else if (fastfabs(p.z())<=fDz+kCarTolerance*0.5)
			{
			    in=kSurface;
			}
		}
	    else if (fastfabs(p.y())<=fDy+kCarTolerance*0.5)
		{
		    if (fastfabs(p.z())<=fDz+kCarTolerance*0.5)
			{
			    in=kSurface;
			}
		}
	}
    else if (fastfabs(p.x())<=fDx+kCarTolerance*0.5)
	{
	    if (fastfabs(p.y())<=fDy+kCarTolerance*0.5)
		{
		    if (fastfabs(p.z())<=fDz+kCarTolerance*0.5)
			{
			    in=kSurface;
			}
		}
	}

    return in;
}

// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned
G4ThreeVector G4Box::SurfaceNormal( const G4ThreeVector& p) const
{
    G4double distx,disty,distz;
    G4ThreeVector norm;
// Calculate distances as if in 1st octant
    distx=fastfabs(fastfabs(p.x())-fDx);
    disty=fastfabs(fastfabs(p.y())-fDy);
    distz=fastfabs(fastfabs(p.z())-fDz);

    if (distx<=disty)
	{
	    if (distx<=distz)
		{
// Closest to X
		    if (p.x()<0) norm=G4ThreeVector(-1.0,0,0);
		    else         norm=G4ThreeVector(1.0,0,0);
		}
	    else
		{
// Closest to Z
		    if (p.z()<0) norm=G4ThreeVector(0,0,-1.0);
		    else         norm=G4ThreeVector(0,0,1.0);
		}
	}
    else
	{
	    if (disty<=distz)
		{
// Closest to Y
		    if (p.y()<0) norm=G4ThreeVector(0,-1.0,0);
		    else         norm=G4ThreeVector(0,1.0,0);
		}
	    else
		{
// Closest to Z
		    if (p.z()<0) norm=G4ThreeVector(0,0,-1.0);
		    else         norm=G4ThreeVector(0,0,1.0);
		}
	}

    return norm;
}

// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possiblity of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is inside the exact
// shape.

G4double G4Box::DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const
{
    G4double safx,safy,safz;
    G4double smin,sminy,sminz;
    G4double smax,smaxy,smaxz;
    G4double stmp;

    safx=fastfabs(p.x())-fDx;   // minimum distance to x surface of shape
    safy=fastfabs(p.y())-fDy;
    safz=fastfabs(p.z())-fDz;

// Will we intersect?
// If safx/y/z is >-tol/2 the point is outside/on the box's x/y/z extent.
// If both p.x/y/z and v.x/y/z repectively are both positive/negative,
// travel is in a direction away from the shape.

    if (   ((p.x()*v.x()>=0.0) && safx>-kCarTolerance*0.5) 
	|| ((p.y()*v.y()>=0.0) && safy>-kCarTolerance*0.5)
        || ((p.z()*v.z()>=0.0) && safz>-kCarTolerance*0.5)) return kInfinity;

// Compute min / max distances for x/y/z travel:
// X Planes
    if (v.x())
        {
            stmp=1.0/fastfabs(v.x());
            smin=safx*stmp;
            smax=(fDx+fastfabs(p.x()))*stmp;
        }
    else
        {
            if (safx<=0.0)
                {
                    smin=0.0;
                    smax=kInfinity;
                }
            else
                {
                    return kInfinity; // Travel parallel
                }
        }

// Y Planes
    if (v.y())
        {
            stmp=1.0/fastfabs(v.y());
            sminy=safy*stmp;
            smaxy=(fDy+fastfabs(p.y()))*stmp;
            if (sminy>smin) smin=sminy;
            if (smaxy<smax) smax=smaxy;
            if (smin>smax) return kInfinity;
        }
    else
        {
            if (safy>0.0)
                {
                    return kInfinity; // Travel parallel
                }
        }

// Z planes
    if (v.z())
        {
            stmp=1.0/fastfabs(v.z());
            sminz=safz*stmp;
            smaxz=(fDz+fastfabs(p.z()))*stmp;
            if (sminz>smin) smin=sminz;
            if (smaxz<smax) smax=smaxz;
            if (smin>smax) return kInfinity;
        }
    else
        {
            if (safz>0.0)
                {
                    return kInfinity; // Travel parallel
                }
        }


    if (smin<0)
        {
            return 0.0;
        }
    return smin;
}
 
// Appoximate distance to box.
// Returns largest perpendicular distance to the closest x/y/z sides of
// the box.
// - If inside return 0
G4double G4Box::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safex,safey,safez,safe=0.0;
    safex=fastfabs(p.x())-fDx;
    safey=fastfabs(p.y())-fDy;
    safez=fastfabs(p.z())-fDz;

    if (safex>safe) safe=safex;
    if (safey>safe) safe=safey;
    if (safez>safe) safe=safez;
    return safe;
}

// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.
// - Eliminate one side of each pair by considering direction of v
// - when leaving a surface & v.close, return 0
G4double G4Box::DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			      const G4bool calcNorm,
			      G4bool *validNorm,G4ThreeVector *n) const
{
    ESide side;
    G4double pdist,stmp,snxt;

    if (calcNorm) *validNorm=true; // All normals are valid

    if (v.x()>0)
	{
	    pdist=fDx-p.x();
	    if (pdist>kCarTolerance*0.5)
		{
		    snxt=pdist/v.x();
		    side=kPX;
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(1,0,0);
			}
		    return snxt=0;
		}
	}
    else if (v.x()<0) 
	{
	    pdist=fDx+p.x();
	    if (pdist>kCarTolerance*0.5)
		{
		    snxt=-pdist/v.x();
		    side=kMX;
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(-1,0,0);
			}
		    return snxt=0;
		}
	}
    else
	{
	    snxt=kInfinity;
	}

    if (v.y()>0)
	{
	    pdist=fDy-p.y();
	    if (pdist>kCarTolerance*0.5)
		{
		    stmp=pdist/v.y();
		    if (stmp<snxt)
			{
			    snxt=stmp;
			    side=kPY;
			}
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(0,1,0);
			}
		    return snxt=0;
		}
	}
    else if (v.y()<0) 
	{
	    pdist=fDy+p.y();
	    if (pdist>kCarTolerance*0.5)
		{
		    stmp=-pdist/v.y();
		    if (stmp<snxt)
			{
			    snxt=stmp;
			    side=kMY;
			}
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(0,-1,0);
			}
		    return snxt=0;
		}
	}

    if (v.z()>0)
	{
	    pdist=fDz-p.z();
	    if (pdist>kCarTolerance*0.5)
		{
		    stmp=pdist/v.z();
		    if (stmp<snxt)
			{
			    snxt=stmp;
			    side=kPZ;
			}
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(0,0,1);
			}
		    return snxt=0;
		}
	}
    else if (v.z()<0) 
	{
	    pdist=fDz+p.z();
	    if (pdist>kCarTolerance*0.5)
		{
		    stmp=-pdist/v.z();
		    if (stmp<snxt)
			{
			    snxt=stmp;
			    side=kMZ;
			}
		}
	    else
		{
		    if (calcNorm)
			{
			    *n=G4ThreeVector(0,0,-1);
			}
		    return snxt=0;
		}
	}


    if (calcNorm)
	{
	    
	    switch (side)
		{
		case kPX:
		    *n=G4ThreeVector(1,0,0);
		    break;
		case kMX:
		    *n=G4ThreeVector(-1,0,0);
		    break;
		case kPY:
		    *n=G4ThreeVector(0,1,0);
		    break;
		case kMY:
		    *n=G4ThreeVector(0,-1,0);
		    break;
		case kPZ:
		    *n=G4ThreeVector(0,0,1);
		    break;
		case kMZ:
		    *n=G4ThreeVector(0,0,-1);
		    break;
		}
	}
    return snxt;
}

// Calculate exact shortest distance to any boundary from inside
// - If outside return 0
G4double G4Box::DistanceToOut(const G4ThreeVector& p) const
{
    G4double safx1,safx2,safy1,safy2,safz1,safz2,safe;
	
    safx1=fDx-p.x();
    safx2=fDx+p.x();
    safy1=fDy-p.y();
    safy2=fDy+p.y();
    safz1=fDz-p.z();
    safz2=fDz+p.z();	
	
// shortest Dist to any boundary now MIN(safx1,safx2,safy1..)
    if (safx2<safx1) safe=safx2;
    else             safe=safx1;
    if (safy1<safe) safe=safy1;
    if (safy2<safe) safe=safy2;
    if (safz1<safe) safe=safz1;
    if (safz2<safe) safe=safz2;

    if (safe<0) safe=0;
    return safe;	
}


