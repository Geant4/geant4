// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trap.cc,v 1.8 2000-11-20 17:58:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Trap
//
// Implementation for G4Trap class
//
// History:
// 21.03.95 P.Kent: Modified for `tolerant' geometry
//  9.09.96 V. Grichine: Final modifications before to commit
//  1.11.96 V.Grichine Costructor for Right Angular Wedge from STEP & G4Trd/Para
//  8.12.97 J.Allison: Added "nominal" constructor and method SetAllParameters.
//  4.06.99 S.Giani: Fixed CalculateExtent in rotated case. 
// 19.11.99 V.Grichine, kUndefined was added to Eside enum
// 13.12.99 V.Grichine, bug fixed in DistanceToIn(p,v)

#include <math.h>
#include "G4Trap.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"

////////////////////////////////////////////////////////////////////////
//
// Accuracy of coplanarity

const G4double kCoplanar_Tolerance=1E-4;

//////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use 
  	
enum Eside {kUndefined,ks0,ks1,ks2,ks3,kPZ,kMZ};

////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Trap::~G4Trap()
{
    ;
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor - check and set half-widths as well as angles: 
// final check of coplanarity

G4Trap::G4Trap( const G4String& pName,
	        G4double pDz,
	        G4double pTheta, G4double pPhi,
	        G4double pDy1, G4double pDx1, G4double pDx2,
		G4double pAlp1,
	        G4double pDy2, G4double pDx3, G4double pDx4,
	        G4double pAlp2) : G4CSGSolid(pName)
{
   if (pDz>0 && pDy1>0 && pDx1>0 && pDx2>0 && pDy2>0 && pDx3>0 && pDx4>0)
   {
	   fDz=pDz;
	   fTthetaCphi=tan(pTheta)*cos(pPhi);
	   fTthetaSphi=tan(pTheta)*sin(pPhi);
	    
	   fDy1=pDy1;
	   fDx1=pDx1;
	   fDx2=pDx2;
	   fTalpha1=tan(pAlp1);
	   
	   fDy2=pDy2;
	   fDx3=pDx3;
	   fDx4=pDx4;
	   fTalpha2=tan(pAlp2);

	   MakePlanes();
    }
    else
    {
       G4Exception("Error in G4Trap::G4Trap - Invalid Length G4Trap parameters");
    }
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor - Design of trapezoid based on 8 G4ThreeVector parameters, 
// which are its vertices. Checking of planarity with preparation of 
// fPlanes[] and than calculation of other members

G4Trap::G4Trap( const G4String& pName,
	        const G4ThreeVector pt[8]) : G4CSGSolid(pName)
{
    if (   pt[0].z()<0 && pt[0].z()==pt[1].z() && pt[0].z()==pt[2].z() && pt[0].z()==pt[3].z()
	&& pt[4].z()>0 && pt[4].z()==pt[5].z() && pt[4].z()==pt[6].z() && pt[4].z()==pt[7].z()
	&& (pt[0].z()+pt[4].z())== 0
	&& pt[0].y()==pt[1].y() && pt[2].y()==pt[3].y()
	&& pt[4].y()==pt[5].y() && pt[6].y()==pt[7].y()
	&& (pt[0].y()+pt[2].y()+pt[4].y()+pt[6].y())==0      )
     {
        G4bool good;
    
// Bottom side with normal approx. -Y
	    good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-Y not planar");

		}
// Top side with normal approx. +Y
	    good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+Y not planar");
		}
// Front side with normal approx. -X
	    good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-X not planar");
		}
// Back side iwth normal approx. +X
	    good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+X not planar");
		}

	   fDz = (pt[7]).z() ;
	    
	   fDy1 = ((pt[2]).y()-(pt[1]).y())*0.5 ;
	   fDx1 = ((pt[1]).x()-(pt[0]).x())*0.5 ;
	   fDx2 = ((pt[3]).x()-(pt[2]).x())*0.5 ;
	   fTalpha1 = ((pt[2]).x()+(pt[3]).x()-(pt[1]).x()-(pt[0]).x())*0.25/fDy1 ;

	   fDy2 = ((pt[6]).y()-(pt[5]).y())*0.5 ;
	   fDx3 = ((pt[5]).x()-(pt[4]).x())*0.5 ;
	   fDx4 = ((pt[7]).x()-(pt[6]).x())*0.5 ;
	   fTalpha2 = ((pt[6]).x()+(pt[7]).x()-(pt[5]).x()-(pt[4]).x())*0.25/fDy2 ;

	   fTthetaCphi = ((pt[4]).x()+fDy2*fTalpha2+fDx3)/fDz ;
	   fTthetaSphi = ((pt[4]).y()+fDy2)/fDz ;

	}
    else
	{
	    G4Exception("Error in G4Trap::G4Trap - Invalid vertice coordinates");
	}

    
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor for Right Angular Wedge from STEP

G4Trap::G4Trap( const G4String& pName,
	        G4double pZ,
	        G4double pY,
		G4double pX, G4double pLTX) : G4CSGSolid(pName) 
	       
{
    G4bool good;

    if (pZ>0 && pY>0 && pX>0 && pLTX>0 && pLTX<=pX)
	{
	   fDz = 0.5*pZ ;
	   fTthetaCphi = 0 ;
	   fTthetaSphi = 0 ;
	    
	   fDy1 = 0.5*pY;
	   fDx1 = 0.5*pX ;
	   fDx2 = 0.5*pLTX;
	   fTalpha1 =  0.5*(pLTX - pX)/pY;
	   
	   fDy2 = fDy1 ;
	   fDx3 = fDx1;
	   fDx4 = fDx2 ;
	   fTalpha2 = fTalpha1 ;

	   G4ThreeVector pt[8] ;
	   
	   pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
           pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	   pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
	   pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
// Bottom side with normal approx. -Y
	    good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-Y not planar");

		}
// Top side with normal approx. +Y
	    good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+Y not planar");
		}
// Front side with normal approx. -X
	    good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-X not planar");
		}
// Back side iwth normal approx. +X
	    good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+X not planar");
		}
	}
    else
	{
	    G4Exception("Error in G4Trap::G4Trap - Invalid Length Trapmeters");
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Trd

G4Trap::G4Trap( const G4String& pName,
	        G4double pDx1,  G4double pDx2,
	        G4double pDy1,  G4double pDy2,
		G4double pDz) : G4CSGSolid(pName)
{
    G4bool good;

    if (pDz>0 && pDy1>0 && pDx1>0 && pDx2>0 && pDy2>0 )
	{
	   fDz = pDz;
	   fTthetaCphi = 0 ;
	   fTthetaSphi = 0 ;
	    
	   fDy1 = pDy1 ;
	   fDx1 = pDx1 ;
	   fDx2 = pDx1 ;
	   fTalpha1 = 0 ;
	   
	   fDy2 = pDy2 ;
	   fDx3 = pDx2 ;
	   fDx4 = pDx2 ;
	   fTalpha2 = 0 ;

	   G4ThreeVector pt[8] ;
	   
	   pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
           pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	   pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
	   pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
// Bottom side with normal approx. -Y
	    good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-Y not planar");

		}
// Top side with normal approx. +Y
	    good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+Y not planar");
		}
// Front side with normal approx. -X
	    good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-X not planar");
		}
// Back side iwth normal approx. +X
	    good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+X not planar");
		}
	}
    else
	{
	    G4Exception("Error in G4Trap::G4Trap - Invalid Length G4Trapmeters");
	}
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor for G4Para

G4Trap::G4Trap( const G4String& pName,
		G4double pDx, G4double pDy,
	        G4double pDz,
	        G4double pAlpha,
		G4double pTheta, G4double pPhi) : G4CSGSolid(pName)
	       
{
    G4bool good;

    if (pDz>0 && pDy>0 && pDx>0 )
	{
	   fDz = pDz ;
	   fTthetaCphi = tan(pTheta)*cos(pPhi) ;
	   fTthetaSphi = tan(pTheta)*sin(pPhi) ;
	    
	   fDy1 = pDy ;
	   fDx1 = pDx ;
	   fDx2 = pDx ;
	   fTalpha1 = tan(pAlpha) ;
	   
	   fDy2 = pDy ;
	   fDx3 = pDx ;
	   fDx4 = pDx ;
	   fTalpha2 = fTalpha1 ;

	   G4ThreeVector pt[8] ;
	   
	   pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
           pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	   pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
	   pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
// Bottom side with normal approx. -Y
	    good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-Y not planar");

		}
// Top side with normal approx. +Y
	    good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+Y not planar");
		}
// Front side with normal approx. -X
	    good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~-X not planar");
		}
// Back side iwth normal approx. +X
	    good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
	    if (!good)
		{
		    G4Exception("G4Trap::G4Trap - face at ~+X not planar");
		}
	}
    else
	{
	    G4Exception("Error in G4Trap::G4Trap - Invalid Length G4Trapmeters");
	}
}

///////////////////////////////////////////////////////////////////////////
//
// Nominal constructor for G4Trap whose parameters are to be set by
// a G4VParamaterisation later.  Check and set half-widths as well as
// angles: final check of coplanarity

G4Trap::G4Trap( const G4String& pName) :
  G4CSGSolid (pName),
  fDz         (1.),
  fTthetaCphi (0.),
  fTthetaSphi (0.),
  fDy1        (1.),
  fDx1        (1.),
  fDx2        (1.),
  fTalpha1    (0.),
  fDy2        (1.),
  fDx3        (1.),
  fDx4        (1.),
  fTalpha2    (0.)
{
  ;
}

///////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - check and set half-widths
// as well as angles: final check of coplanarity

void G4Trap::SetAllParameters ( G4double pDz,
				G4double pTheta,
				G4double pPhi,
				G4double pDy1,
				G4double pDx1,
				G4double pDx2,
				G4double pAlp1,
				G4double pDy2,
				G4double pDx3,
				G4double pDx4,
				G4double pAlp2)
{
    if (pDz>0 && pDy1>0 && pDx1>0 && pDx2>0 && pDy2>0 && pDx3>0 && pDx4>0)
	{
	   fDz=pDz;
	   fTthetaCphi=tan(pTheta)*cos(pPhi);
	   fTthetaSphi=tan(pTheta)*sin(pPhi);
	    
	   fDy1=pDy1;
	   fDx1=pDx1;
	   fDx2=pDx2;
	   fTalpha1=tan(pAlp1);
	   
	   fDy2=pDy2;
	   fDx3=pDx3;
	   fDx4=pDx4;
	   fTalpha2=tan(pAlp2);

	   MakePlanes();
	}
    else
    {
       G4Exception("Error in G4Trap::SetAll - Invalid Length G4Trap parameters");
    }
}

//////////////////////////////////////////////////////////////////////////
//
// Checking of coplanarity

G4bool G4Trap::MakePlanes()
{
  G4bool good = true;

  G4ThreeVector pt[8] ;
	   
  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
  // Bottom side with normal approx. -Y

  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]) ;

  if (!good) G4Exception("G4Trap::G4Trap - face at ~-Y not planar");

  // Top side with normal approx. +Y

  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);

  if (!good) G4Exception("G4Trap::G4Trap - face at ~+Y not planar");

  // Front side with normal approx. -X

  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);

  if (!good) G4Exception("G4Trap::G4Trap - face at ~-X not planar");
   
  // Back side iwth normal approx. +X

  good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);

  if (!good) G4Exception("G4Trap::G4Trap - face at ~+X not planar");
    
  return good;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate the coef's of the plane p1->p2->p3->p4->p1
// where the ThreeVectors 1-4 are in anti-clockwise order when viewed from
// infront of the plane.
//
// Return true if the ThreeVectors are coplanar + set coef;s
//        false if ThreeVectors are not coplanar

G4bool G4Trap::MakePlane( const G4ThreeVector& p1,
			  const G4ThreeVector& p2,
		          const G4ThreeVector& p3,
			  const G4ThreeVector& p4,
			  TrapSidePlane& plane )
{
    G4double a,b,c,s;
    G4ThreeVector v12,v13,v14,Vcross;

    G4bool good;

    v12 = p2-p1;
    v13 = p3-p1;
    v14 = p4-p1;
    Vcross=v12.cross(v13);

    if (fabs(Vcross.dot(v14)/(Vcross.mag()*v14.mag())) > kCoplanar_Tolerance)
    {
       good=false;
    }
    else
    {
    
// a,b,c correspond to the x/y/z components of the normal vector to the plane
	   
            a=(p2.y()-p1.y())*(p1.z()+p2.z())+(p3.y()-p2.y())*(p2.z()+p3.z());
	    a+=(p4.y()-p3.y())*(p3.z()+p4.z())+(p1.y()-p4.y())*(p4.z()+p1.z()); // may be delete ?
    
	    b=(p2.z()-p1.z())*(p1.x()+p2.x())+(p3.z()-p2.z())*(p2.x()+p3.x());
	    b+=(p4.z()-p3.z())*(p3.x()+p4.x())+(p1.z()-p4.z())*(p4.x()+p1.x()); // ?
	    
	    c=(p2.x()-p1.x())*(p1.y()+p2.y())+(p3.x()-p2.x())*(p2.y()+p3.y());
	    c+=(p4.x()-p3.x())*(p3.y()+p4.y())+(p1.x()-p4.x())*(p4.y()+p1.y());  // ?
// Let create diagonals 4-2 and 3-1 than (4-2)x(3-1) provides vector perpendicular to the
// plane directed to outside !!! and a,b,c, = f(1,2,3,4)	   
	   //a = +(p4.y() - p2.y())*(p3.z() - p1.z()) - (p3.y() - p1.y())*(p4.z() - p2.z()) ;
	   //b = -(p4.x() - p2.x())*(p3.z() - p1.z()) + (p3.x() - p1.x())*(p4.z() - p2.z()) ; 
	   //c = +(p4.x() - p2.x())*(p3.y() - p1.y()) - (p3.x() - p1.x())*(p4.y() - p2.y()) ;
	   s=sqrt(a*a+b*b+c*c);   // so now vector plane.(a,b,c) is unit 

	   if( s > 0 )
	   {
		plane.a=a/s;
		plane.b=b/s;
		plane.c=c/s;
	   }
	   else 
                G4Exception("Invalid parameters in G4Trap::MakePlane") ;

// Calculate D: p1 in in plane so D=-n.p1.Vect()
	   
   plane.d=-(plane.a*p1.x()+plane.b*p1.y()+plane.c*p1.z());
   good=true;
 }
 return good;
}

//////////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

 void G4Trap::ComputeDimensions(G4VPVParameterisation* p,
                                const G4int n,
                                const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}


////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Trap::CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                               G4double& pMin, G4double& pMax) const
{

    G4double xMin, xMax, yMin, yMax, zMin, zMax;
    G4bool flag;

    if (!pTransform.IsRotated())
   {  
                      // Special case handling for unrotated trapezoids
                      // Compute z/x/y/ mins and maxs respecting limits, with early returns
                      // if outside limits. Then switch() on pAxis
	  G4int i ; 
          G4double xoffset;
          G4double yoffset;
          G4double zoffset;
	  G4double temp[8] ;          // some points for intersection with zMin/zMax
	  
          xoffset=pTransform.NetTranslation().x();	    
          yoffset=pTransform.NetTranslation().y();
          zoffset=pTransform.NetTranslation().z();
 
          G4ThreeVector pt[8];   // vertices after translation
          pt[0]=G4ThreeVector(xoffset-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
               		      yoffset-fDz*fTthetaSphi-fDy1,zoffset-fDz);
          pt[1]=G4ThreeVector(xoffset-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
		              yoffset-fDz*fTthetaSphi-fDy1,zoffset-fDz);
          pt[2]=G4ThreeVector(xoffset-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
		              yoffset-fDz*fTthetaSphi+fDy1,zoffset-fDz);
          pt[3]=G4ThreeVector(xoffset-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
		              yoffset-fDz*fTthetaSphi+fDy1,zoffset-fDz);
          pt[4]=G4ThreeVector(xoffset+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
		              yoffset+fDz*fTthetaSphi-fDy2,zoffset+fDz);
          pt[5]=G4ThreeVector(xoffset+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
		              yoffset+fDz*fTthetaSphi-fDy2,zoffset+fDz);
          pt[6]=G4ThreeVector(xoffset+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
		              yoffset+fDz*fTthetaSphi+fDy2,zoffset+fDz);
          pt[7]=G4ThreeVector(xoffset+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
		              yoffset+fDz*fTthetaSphi+fDy2,zoffset+fDz);
	    zMin=zoffset-fDz;
	    zMax=zoffset+fDz;
	    if (pVoxelLimit.IsZLimited())
		{
		    if (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance
			||zMax<pVoxelLimit.GetMinZExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (zMin<pVoxelLimit.GetMinZExtent())
				{
				    zMin=pVoxelLimit.GetMinZExtent();
				}
			    if (zMax>pVoxelLimit.GetMaxZExtent())
				{
				    zMax=pVoxelLimit.GetMaxZExtent();
				}
			}
		}

            temp[0] = pt[0].y()+(pt[4].y()-pt[0].y())*(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
       	    temp[1] = pt[0].y()+(pt[4].y()-pt[0].y())*(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
	    temp[2] = pt[2].y()+(pt[6].y()-pt[2].y())*(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
	    temp[3] = pt[2].y()+(pt[6].y()-pt[2].y())*(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;	      
	    yMax = yoffset - fabs(fDz*fTthetaSphi) - fDy1 - fDy2 ;
	    yMin = -yMax ;
	    for(i=0;i<4;i++)
	    {
	       if(temp[i] > yMax) yMax = temp[i] ;
	       if(temp[i] < yMin) yMin = temp[i] ;
	    }
	    
	    if (pVoxelLimit.IsYLimited())
		{
		    if (yMin>pVoxelLimit.GetMaxYExtent()+kCarTolerance
			||yMax<pVoxelLimit.GetMinYExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (yMin<pVoxelLimit.GetMinYExtent())
				{
				    yMin=pVoxelLimit.GetMinYExtent();
				}
			    if (yMax>pVoxelLimit.GetMaxYExtent())
				{
				    yMax=pVoxelLimit.GetMaxYExtent();
				}
			}
		}

            temp[0] = pt[0].x()+(pt[4].x()-pt[0].x())*(zMin-pt[0].z())/(pt[4].z()-pt[0].z()) ;
       	    temp[1] = pt[0].x()+(pt[4].x()-pt[0].x())*(zMax-pt[0].z())/(pt[4].z()-pt[0].z()) ;
	    temp[2] = pt[2].x()+(pt[6].x()-pt[2].x())*(zMin-pt[2].z())/(pt[6].z()-pt[2].z()) ;
	    temp[3] = pt[2].x()+(pt[6].x()-pt[2].x())*(zMax-pt[2].z())/(pt[6].z()-pt[2].z()) ;
            temp[4] = pt[3].x()+(pt[7].x()-pt[3].x())*(zMin-pt[3].z())/(pt[7].z()-pt[3].z()) ;
       	    temp[5] = pt[3].x()+(pt[7].x()-pt[3].x())*(zMax-pt[3].z())/(pt[7].z()-pt[3].z()) ;
	    temp[6] = pt[1].x()+(pt[5].x()-pt[1].x())*(zMin-pt[1].z())/(pt[5].z()-pt[1].z()) ;
	    temp[7] = pt[1].x()+(pt[5].x()-pt[1].x())*(zMax-pt[1].z())/(pt[5].z()-pt[1].z()) ;
	    
	    xMax = xoffset - fabs(fDz*fTthetaCphi) - fDx1 - fDx2 -fDx3 - fDx4 ;
	    xMin = -xMax ;
	    for(i=0;i<8;i++)
	    {
	       if(temp[i] > xMax) xMax = temp[i] ;
	       if(temp[i] < xMin) xMin = temp[i] ;
	    }
	                                          // xMax/Min = f(yMax/Min) ?
	    if (pVoxelLimit.IsXLimited())
		{
		    if (xMin>pVoxelLimit.GetMaxXExtent()+kCarTolerance
			||xMax<pVoxelLimit.GetMinXExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (xMin<pVoxelLimit.GetMinXExtent())
				{
				    xMin=pVoxelLimit.GetMinXExtent();
				}
			    if (xMax>pVoxelLimit.GetMaxXExtent())
				{
				    xMax=pVoxelLimit.GetMaxXExtent();
				}
			}
		}

	    switch (pAxis)
		{
		case kXAxis:
		    pMin=xMin;
		    pMax=xMax;
		    break;
		case kYAxis:
		    pMin=yMin;
		    pMax=yMax;
		    break;
		case kZAxis:
		    pMin=zMin;
		    pMax=zMax;
		    break;
		default:
		    break;
		}

	    pMin-=kCarTolerance;
	    pMax+=kCarTolerance;

	    flag = true;
   }
    else
   {
// General rotated case - 

	    G4bool existsAfterClip=false;
	    G4ThreeVectorList *vertices;

	    pMin=+kInfinity;
	    pMax=-kInfinity;
	    
// Calculate rotated vertex coordinates
	    vertices=CreateRotatedVertices(pTransform);
	    
	    xMin = +kInfinity; yMin = +kInfinity; zMin = +kInfinity;
	    xMax = -kInfinity; yMax = -kInfinity; zMax = -kInfinity;
	    
	    for(G4int nv=0; nv<8; nv++){
	    
	       if((*vertices)[nv].x() > xMax){xMax = (*vertices)[nv].x();};
	       if((*vertices)[nv].y() > yMax){yMax = (*vertices)[nv].y();};
	       if((*vertices)[nv].z() > zMax){zMax = (*vertices)[nv].z();};
	    
	       if((*vertices)[nv].x() < xMin){xMin = (*vertices)[nv].x();};
	       if((*vertices)[nv].y() < yMin){yMin = (*vertices)[nv].y();};
	       if((*vertices)[nv].z() < zMin){zMin = (*vertices)[nv].z();};
	       
	    };
	    
	    if (pVoxelLimit.IsZLimited())
		{
		    if (zMin>pVoxelLimit.GetMaxZExtent()+kCarTolerance
			||zMax<pVoxelLimit.GetMinZExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (zMin<pVoxelLimit.GetMinZExtent())
				{
				    zMin=pVoxelLimit.GetMinZExtent();
				}
			    if (zMax>pVoxelLimit.GetMaxZExtent())
				{
				    zMax=pVoxelLimit.GetMaxZExtent();
				}
			}
		}
	    
            if (pVoxelLimit.IsYLimited())
		{
		    if (yMin>pVoxelLimit.GetMaxYExtent()+kCarTolerance
			||yMax<pVoxelLimit.GetMinYExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (yMin<pVoxelLimit.GetMinYExtent())
				{
				    yMin=pVoxelLimit.GetMinYExtent();
				}
			    if (yMax>pVoxelLimit.GetMaxYExtent())
				{
				    yMax=pVoxelLimit.GetMaxYExtent();
				}
			}
		}

             if (pVoxelLimit.IsXLimited())
		{
		    if (xMin>pVoxelLimit.GetMaxXExtent()+kCarTolerance
			||xMax<pVoxelLimit.GetMinXExtent()-kCarTolerance)
			{
			    return false;
			}
		    else
			{
			    if (xMin<pVoxelLimit.GetMinXExtent())
				{
				    xMin=pVoxelLimit.GetMinXExtent();
				}
			    if (xMax>pVoxelLimit.GetMaxXExtent())
				{
				    xMax=pVoxelLimit.GetMaxXExtent();
				}
			}
		}

	    switch (pAxis)
		{
		case kXAxis:
		    pMin=xMin;
		    pMax=xMax;
		    break;
		case kYAxis:
		    pMin=yMin;
		    pMax=yMax;
		    break;
		case kZAxis:
		    pMin=zMin;
		    pMax=zMax;
		    break;
		default:
		    break;
		}

	    
	    if (pMin!=kInfinity||pMax!=-kInfinity)
		{
		    existsAfterClip=true;
		    
// Add tolerance to avoid precision troubles
		    pMin-=kCarTolerance;
		    pMax+=kCarTolerance;
		    
		};

	    delete vertices ;          //  'new' in the function called
	    flag = existsAfterClip ;
   }
   return flag;
}


////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Trap::Inside(const G4ThreeVector& p) const
{
    EInside in;
    G4double Dist;
    G4int i;
    if (fabs(p.z())<=fDz-kCarTolerance/2)
	{
	    in=kInside;
	    for (i=0;i<4;i++)
		{
		    Dist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d;
		    if (Dist>kCarTolerance/2)
			{
			    return in=kOutside;
			}
		    else if (Dist>-kCarTolerance/2)
			{
			    in=kSurface;
			} 
		}
	}
    else if (fabs(p.z())<=fDz+kCarTolerance/2)
	{
	    in=kSurface;
	    for (i=0;i<4;i++)
		{
		    Dist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d;
		    if (Dist>kCarTolerance/2)
			{
			    return in=kOutside;
			}
		}

	}
    else
	{
	    in=kOutside;
	}
    return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)

G4ThreeVector G4Trap::SurfaceNormal( const G4ThreeVector& p) const
{
    G4double safe=kInfinity,Dist,safez;
    G4int i,imin=0;
    for (i=0;i<4;i++)
	{
	    Dist=fabs(fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d);
	    if (Dist<safe)
		{
		   safe=Dist;
		    imin=i;
		}
	}

    safez=fabs(fabs(p.z())-fDz);
    if (safe<safez)
	{
	    return G4ThreeVector(fPlanes[imin].a,fPlanes[imin].b,fPlanes[imin].c);
	}
    else
	{
	    if (p.z()>0)
		{
		    return G4ThreeVector(0,0,1);
		}
	    else
		{
		    return G4ThreeVector(0,0,-1);
		}
	}
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside - return kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect

G4double G4Trap::DistanceToIn(const G4ThreeVector& p,
			      const G4ThreeVector& v) const
{

    G4double snxt;		// snxt = default return value
    G4double max,smax,smin;
    G4double pdist,Comp,vdist;
    G4int i;
//
// Z Intersection range
//
    if (v.z() > 0 )
    {
       max = fDz - p.z() ;

       if (max > 0.5*kCarTolerance)
       {
	  smax = max/v.z();
	  smin = (-fDz-p.z())/v.z();
       }
       else
       {
	  return snxt=kInfinity;
       }
    }
    else if (v.z() < 0 )
    {
       max = - fDz - p.z() ;

       if (max < -0.5*kCarTolerance )
       {
	  smax=max/v.z();
	  smin=(fDz-p.z())/v.z();
       }
       else
       {
	  return snxt=kInfinity;
       }
    }
    else
    {
       if (fabs(p.z())<fDz - 0.5*kCarTolerance) // Inside was <=fDz
       {
	  smin=0;
	  smax=kInfinity;
       }
       else
       {
	  return snxt=kInfinity;
       }
    }

    for (i=0;i<4;i++)
    {
       pdist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d;
       
       Comp=fPlanes[i].a*v.x()+fPlanes[i].b*v.y()+fPlanes[i].c*v.z();

       if ( pdist >= -0.5*kCarTolerance )      // was >0
       {
//
// Outside the plane -> this is an extent entry distance
//
	  if (Comp >= 0)   // was >0
	  {
	     return snxt=kInfinity ;
	  }
	  else 
	  {
	     vdist=-pdist/Comp;

	     if (vdist>smin)
	     {
		if (vdist<smax)
		{
		   smin = vdist;
		}
		else
		{
		   return snxt=kInfinity;
		}
	     }
	  }
       }
       else
       {
//
// Inside the plane -> couble  be an extent exit distance (smax)
//
	  if (Comp>0)  // Will leave extent
	  {
	     vdist=-pdist/Comp;

	     if (vdist<smax)
	     {
		if (vdist>smin)
		{
		   smax=vdist;
		}
		else
		{
		       return snxt=kInfinity;
		}
	     }	
	  }
       }
   }
//
// Checks in non z plane intersections ensure smin<smax
//
    if (smin >=0 )
    {
       snxt = smin ;
    }
    else
    {
       snxt = 0 ;
    }
    return snxt;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - Returns 0 is ThreeVector inside

G4double G4Trap::DistanceToIn(const G4ThreeVector& p) const
{
    G4double safe,Dist;
    G4int i;
    safe=fabs(p.z())-fDz;
    for (i=0;i<4;i++)
    {
      Dist=fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d;

      if (Dist > safe) safe=Dist;
    }
    if (safe<0) safe=0;
    return safe;	
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance

G4double G4Trap::DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			       const G4bool calcNorm,
			       G4bool *validNorm,G4ThreeVector *n) const
{
    Eside side = kUndefined ;
    G4double snxt;		// snxt = return value
    G4double pdist,Comp,vdist,max;
//
// Z Intersections
//
    if (v.z()>0)
	{
	    max=fDz-p.z();
	    if (max>kCarTolerance/2)
		{
		   snxt=max/v.z();
		   side=kPZ;
		}
	    else
		{
		    if (calcNorm)
			{
			    *validNorm=true;
			    *n=G4ThreeVector(0,0,1);
			}
		    return snxt=0;
		}
	}
    else if (v.z()<0)
	{
	    max=-fDz-p.z();
	    if (max<-kCarTolerance/2)
		{
		   snxt=max/v.z();
		   side=kMZ;
		}
	    else
		{
		    if (calcNorm)
			{
			    *validNorm=true;
			    *n=G4ThreeVector(0,0,-1);
			}
		    return snxt=0;
		}
	}
    else
    {
       snxt=kInfinity;
    }

//
// Intersections with planes[0] (expanded because of setting enum)
//
    pdist=fPlanes[0].a*p.x()+fPlanes[0].b*p.y()+fPlanes[0].c*p.z()+fPlanes[0].d;
    Comp=fPlanes[0].a*v.x()+fPlanes[0].b*v.y()+fPlanes[0].c*v.z();
    if (pdist>0)
	{
// Outside the plane
	    if (Comp>0)
		{
// Leaving immediately
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
		    }
		    return snxt=0;
		}
	}
    else if (pdist<-kCarTolerance/2)
	{
// Inside the plane
	    if (Comp>0)
		{
// Will leave extent
		    vdist=-pdist/Comp;
		    if (vdist<snxt)
			{
			   snxt=vdist;
			   side=ks0;
			}
		}
	}
    else
	{
// On surface
	    if (Comp>0)
		{
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
		    }
		    return snxt=0;
		}
	}


//
// Intersections with planes[1] (expanded because of setting enum)
//
    pdist=fPlanes[1].a*p.x()+fPlanes[1].b*p.y()+fPlanes[1].c*p.z()+fPlanes[1].d;
    Comp=fPlanes[1].a*v.x()+fPlanes[1].b*v.y()+fPlanes[1].c*v.z();
    if (pdist>0)
	{
// Outside the plane
	    if (Comp>0)
		{
// Leaving immediately
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
		    }
		    return snxt=0;
		}
	}
    else if (pdist<-kCarTolerance/2)
	{
// Inside the plane
	    if (Comp>0)
		{
// Will leave extent
		    vdist=-pdist/Comp;
		    if (vdist<snxt)
			{
			   snxt=vdist;
			   side=ks1;
			}
		}
	}
    else
	{
// On surface
	    if (Comp>0)
		{
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
		    }
		    return snxt=0;
		}
	}

//
// Intersections with planes[2] (expanded because of setting enum)
//
    pdist=fPlanes[2].a*p.x()+fPlanes[2].b*p.y()+fPlanes[2].c*p.z()+fPlanes[2].d;
    Comp=fPlanes[2].a*v.x()+fPlanes[2].b*v.y()+fPlanes[2].c*v.z();
    if (pdist>0)
	{
// Outside the plane
	    if (Comp>0)
		{
// Leaving immediately
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
		    }
		    return snxt=0;
		}
	}
    else if (pdist<-kCarTolerance/2)
	{
// Inside the plane
	    if (Comp>0)
		{
// Will leave extent
		    vdist=-pdist/Comp;
		    if (vdist<snxt)
			{
			   snxt=vdist;
			   side=ks2;
			}
		}
	}
    else
	{
// On surface
	    if (Comp>0)
		{
		    if (calcNorm)
		    {
			    *validNorm=true;
		       *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
		    }
		    return snxt=0;
		}
	}

//
// Intersections with planes[3] (expanded because of setting enum)
//
    pdist=fPlanes[3].a*p.x()+fPlanes[3].b*p.y()+fPlanes[3].c*p.z()+fPlanes[3].d;
    Comp=fPlanes[3].a*v.x()+fPlanes[3].b*v.y()+fPlanes[3].c*v.z();
    if (pdist>0)
	{
// Outside the plane
	    if (Comp>0)
		{
// Leaving immediately
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
		    }
		    return snxt=0;
		}
	}
    else if (pdist<-kCarTolerance/2)
	{
// Inside the plane
	    if (Comp>0)
		{
// Will leave extent
		    vdist=-pdist/Comp;
		    if (vdist<snxt)
			{
			   snxt=vdist;
			   side=ks3;
			}
		}
	}
    else
	{
// On surface
	    if (Comp>0)
		{
		    if (calcNorm)
		    {
			    *validNorm=true;
		      *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
		    }
		    return snxt=0;
		}
	}

// set normal
    if (calcNorm)
	{
	    *validNorm=true;
	    switch(side)
	    {
		case ks0:
		    *n=G4ThreeVector(fPlanes[0].a,fPlanes[0].b,fPlanes[0].c);
		    break;
		case ks1:
		    *n=G4ThreeVector(fPlanes[1].a,fPlanes[1].b,fPlanes[1].c);
		    break;
		case ks2:
		    *n=G4ThreeVector(fPlanes[2].a,fPlanes[2].b,fPlanes[2].c);
		    break;
		case ks3:
		    *n=G4ThreeVector(fPlanes[3].a,fPlanes[3].b,fPlanes[3].c);
		    break;
		case kMZ:
		    *n=G4ThreeVector(0,0,-1);
		    break;
		case kPZ:
		    *n=G4ThreeVector(0,0,1);
		    break;
		default:
		    break;
	    }
	}
    return snxt;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is ThreeVector outside

G4double G4Trap::DistanceToOut(const G4ThreeVector& p) const
{
  G4double safe,Dist;
  G4int i;
  safe=fDz-fabs(p.z());

  if (safe<0) safe=0;
  else
  {
    for (i=0;i<4;i++)
    {
     Dist=-(fPlanes[i].a*p.x()+fPlanes[i].b*p.y()+fPlanes[i].c*p.z()+fPlanes[i].d);

      if (Dist<safe) safe=Dist;
    }
    if (safe<0) safe=0;
  }
  return safe;	
}

//////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Ordering [0-3] -fDz cross section
//          [4-7] +fDz cross section such that [0] is below [4],
//                                             [1] below [5] etc.
// Note:
//  Caller has deletion resposibility

G4ThreeVectorList*
G4Trap::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
    G4ThreeVectorList *vertices;
    vertices=new G4ThreeVectorList(8);
    if (vertices)
	{
	    G4ThreeVector vertex0(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	    G4ThreeVector vertex1(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	    G4ThreeVector vertex2(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	    G4ThreeVector vertex3(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	    G4ThreeVector vertex4(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	    G4ThreeVector vertex5(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	    G4ThreeVector vertex6(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
	    G4ThreeVector vertex7(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);

	    vertices->insert(pTransform.TransformPoint(vertex0));
	    vertices->insert(pTransform.TransformPoint(vertex1));
	    vertices->insert(pTransform.TransformPoint(vertex2));
	    vertices->insert(pTransform.TransformPoint(vertex3));
	    vertices->insert(pTransform.TransformPoint(vertex4));
	    vertices->insert(pTransform.TransformPoint(vertex5));
	    vertices->insert(pTransform.TransformPoint(vertex6));
	    vertices->insert(pTransform.TransformPoint(vertex7));
	}
    else
	{
	    G4Exception("G4Trap::CreateRotatedVertices Out of memory - Cannot alloc vertices");
	}
    return vertices;
}

///////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Trap::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddThis (*this);
}

G4Polyhedron* G4Trap::CreatePolyhedron () const
{
    G4double phi = atan2(fTthetaSphi, fTthetaCphi);
    G4double alpha1 = atan(fTalpha1);
    G4double alpha2 = atan(fTalpha2);
    G4double theta = atan(sqrt
                          (fTthetaCphi*fTthetaCphi+fTthetaSphi*fTthetaSphi));
    return new G4PolyhedronTrap(fDz, theta, phi, fDy1, fDx1, fDx2, alpha1,
                                fDy2, fDx3, fDx4, alpha2);
}

G4NURBS* G4Trap::CreateNURBS () const
{
   // return new G4NURBSbox (fDx, fDy, fDz);
   return 0 ;
}


// ********************************  End of G4Trap.cc   ********************************





